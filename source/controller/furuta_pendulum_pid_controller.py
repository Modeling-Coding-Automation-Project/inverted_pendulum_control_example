from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[2]))

import math

TS_DEFAULT = 0.005
TS_MIN = 1e-6

V_LIMIT_DEFAULT = 12.0
ALPHA_REF_LIMIT_RAD_DEFAULT = 12.0 * math.pi / 180.0
THETA_TO_ALPHA_SIGN_DEFAULT = 1.0
KP_THETA_DEFAULT = 0.1
KI_THETA_DEFAULT = 0.0
KD_THETA_DEFAULT = 0.05
KP_ALPHA_DEFAULT = 100.0
KI_ALPHA_DEFAULT = 0.0
KD_ALPHA_DEFAULT = 10.0
DTHETA_TAU_DEFAULT = 0.02
DALPHA_TAU_DEFAULT = 0.02
EPS_VALUE = 1e-12


def _wrap_to_pi(angle_rad: float) -> float:
    """Wrap angle to [-pi, pi]."""
    return (angle_rad + math.pi) % (2.0 * math.pi) - math.pi


class FurutaPendulum_PID_Controller:
    """PID-based controller for a Furuta (rotary inverted) pendulum.

    Inputs:
      - theta, dtheta: arm angle/velocity [rad], [rad/s]
      - alpha, dalpha: pendulum angle/velocity [rad], [rad/s] (alpha=0 is upright)

    References:
      - alpha_ref = 0 rad (upright)
      - theta_ref is stored as a variable so it can be changed later.

    Output:
      - motor voltage command [V]
    """

    def __init__(
        self,
        Ts: float = TS_DEFAULT,
        theta_ref_rad: float = 0.0,
        alpha_ref_rad: float = 0.0,
        v_limit: float = V_LIMIT_DEFAULT,
        alpha_ref_limit_rad: float | None = None,
        theta_to_alpha_sign: float = THETA_TO_ALPHA_SIGN_DEFAULT,
    ) -> None:
        self.Ts = float(Ts)

        # References (can be changed later)
        self.theta_ref_rad = float(theta_ref_rad)
        self.alpha_ref_rad = float(alpha_ref_rad)

        # Voltage saturation used for anti-windup
        self.v_limit = float(abs(v_limit))

        # Limit for commanded pendulum tilt generated from theta control.
        # Keeping this small is important: this controller intentionally uses
        # a *small* pendulum lean to create restoring motion.
        if alpha_ref_limit_rad is None:
            alpha_ref_limit_rad = ALPHA_REF_LIMIT_RAD_DEFAULT
        self.alpha_ref_limit_rad = float(abs(alpha_ref_limit_rad))

        # Mapping sign: if theta is positive, we want alpha_ref_cmd to be positive.
        # If your plant uses opposite sign conventions, set this to -1.0.
        self.theta_to_alpha_sign = 1.0 if theta_to_alpha_sign >= 0.0 else -1.0

        # === Tuning gains (start conservative; tune as needed) ===
        # Theta PID output is interpreted as a commanded pendulum angle offset [rad].
        # (theta > 0) -> (alpha_ref_cmd > 0) by construction.
        self.kp_theta = KP_THETA_DEFAULT
        self.ki_theta = KI_THETA_DEFAULT
        self.kd_theta = KD_THETA_DEFAULT

        # Pendulum PID outputs motor voltage [V] to track alpha_ref_cmd.
        self.kp_alpha = KP_ALPHA_DEFAULT
        self.ki_alpha = KI_ALPHA_DEFAULT
        self.kd_alpha = KD_ALPHA_DEFAULT

        # Integrator states
        self._int_theta = 0.0
        self._int_alpha = 0.0

        # Simple derivative filtering (on measurement)
        self._dalpha_filt = 0.0
        self._dalpha_tau = DALPHA_TAU_DEFAULT  # [s]
        self._dtheta_filt = 0.0
        self._dtheta_tau = DTHETA_TAU_DEFAULT  # [s]

    def reset(self) -> None:
        self._int_theta = 0.0
        self._int_alpha = 0.0
        self._dalpha_filt = 0.0
        self._dtheta_filt = 0.0

    def set_theta_reference_rad(self, theta_ref_rad: float) -> None:
        self.theta_ref_rad = float(theta_ref_rad)

    def set_theta_reference_deg(self, theta_ref_deg: float) -> None:
        self.theta_ref_rad = math.radians(float(theta_ref_deg))

    def set_alpha_reference_rad(self, alpha_ref_rad: float) -> None:
        self.alpha_ref_rad = float(alpha_ref_rad)

    def set_alpha_reference_deg(self, alpha_ref_deg: float) -> None:
        self.alpha_ref_rad = math.radians(float(alpha_ref_deg))

    def calculate_manipulation(
        self,
        theta: float,
        alpha: float,
        dtheta: float,
        dalpha: float,
    ) -> float:
        """Compute voltage command.

                Policy:
                    1) Compute theta PID "manipulated value" (here: desired alpha offset).
                    2) Use that value as the pendulum angle command alpha_ref_cmd.
                    3) Pendulum PID tracks alpha_ref_cmd and outputs motor voltage.

                Intuition:
                    - When theta is positive, alpha_ref_cmd becomes positive (small tilt).
                    - Maintaining that tilt produces a restoring motion on theta through
                        the Furuta coupling, driving both angles back toward 0.
        """

        Ts = self.Ts
        if Ts <= TS_MIN:
            Ts = TS_MIN

        # --- Theta error (defined so theta>ref -> positive error) ---
        theta_meas = float(theta)
        alpha_meas = float(alpha)
        dtheta_meas = float(dtheta)
        dalpha_meas = float(dalpha)

        e_theta = _wrap_to_pi(theta_meas - self.theta_ref_rad)

        # --- Filtered derivatives (measurement) ---
        # Low-pass filter: y_dot_filt = (tau/(tau+Ts))*y_dot_filt + (Ts/(tau+Ts))*y_dot
        tau_th = self._dtheta_tau
        if tau_th > 0.0:
            a = tau_th / (tau_th + Ts)
            b = Ts / (tau_th + Ts)
            self._dtheta_filt = a * self._dtheta_filt + b * dtheta_meas
            dtheta_used = self._dtheta_filt
        else:
            dtheta_used = dtheta_meas

        tau_al = self._dalpha_tau
        if tau_al > 0.0:
            a = tau_al / (tau_al + Ts)
            b = Ts / (tau_al + Ts)
            self._dalpha_filt = a * self._dalpha_filt + b * dalpha_meas
            dalpha_used = self._dalpha_filt
        else:
            dalpha_used = dalpha_meas

        # --- Theta PID -> commanded pendulum angle (alpha_ref_cmd) ---
        # Derivative term uses +dtheta (measurement) so theta increasing positive
        # will also push alpha_ref_cmd positive (consistent with the policy).
        alpha_offset_cmd = self.theta_to_alpha_sign * (
            self.kp_theta * e_theta
            + self.ki_theta * self._int_theta
            + self.kd_theta * dtheta_used
        )

        alpha_ref_cmd = _wrap_to_pi(self.alpha_ref_rad + alpha_offset_cmd)
        alpha_ref_cmd_sat = max(
            -self.alpha_ref_limit_rad,
            min(self.alpha_ref_limit_rad, alpha_ref_cmd),
        )

        # Pendulum error relative to *commanded* reference
        e_alpha = _wrap_to_pi(alpha_ref_cmd_sat - alpha_meas)

        # --- Pendulum stabilization voltage contribution ---
        # Derivative term uses -dalpha (derivative of error when alpha_ref is constant)
        v_unsat = float(
            (self.kp_alpha * e_alpha)
            + (self.ki_alpha * self._int_alpha)
            - (self.kd_alpha * dalpha_used)
        )

        # --- Saturation ---
        v_sat = max(-self.v_limit, min(self.v_limit, v_unsat))

        # --- Anti-windup (conditional integration) ---
        def _should_integrate_u(u_cmd: float, u_cmd_sat: float, err: float, u_lim: float) -> bool:
            if u_cmd == u_cmd_sat:
                return True
            if u_cmd_sat >= u_lim - EPS_VALUE and err > 0.0:
                return False
            if u_cmd_sat <= -u_lim + EPS_VALUE and err < 0.0:
                return False
            return True

        # Theta integrator is limited by the alpha_ref_cmd saturation.
        if _should_integrate_u(alpha_ref_cmd, alpha_ref_cmd_sat, e_theta, self.alpha_ref_limit_rad):
            self._int_theta += e_theta * Ts

        # Alpha integrator is limited by the voltage saturation.
        if _should_integrate_u(v_unsat, v_sat, e_alpha, self.v_limit):
            self._int_alpha += e_alpha * Ts

        return v_sat
