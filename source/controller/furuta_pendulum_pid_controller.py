from __future__ import annotations

import math


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
        Ts: float = 0.005,
        theta_ref_rad: float = 0.0,
        alpha_ref_rad: float = 0.0,
        v_limit: float = 12.0,
    ) -> None:
        self.Ts = float(Ts)

        # References (can be changed later)
        self.theta_ref_rad = float(theta_ref_rad)
        self.alpha_ref_rad = float(alpha_ref_rad)

        # Voltage saturation used for anti-windup
        self.v_limit = float(abs(v_limit))

        # === Tuning gains (start conservative; tune as needed) ===
        # Pendulum (stabilization) PID
        self.kp_alpha = 35.0
        self.ki_alpha = 0.0
        self.kd_alpha = 2.5

        # Arm angle cascade: outer loop (theta -> dtheta_ref)
        self.kp_theta_outer = 3.0
        self.ki_theta_outer = 0.2

        # Inner loop (dtheta -> voltage)
        self.kp_dtheta_inner = 0.25
        self.ki_dtheta_inner = 0.0

        # Integrator states
        self._int_alpha = 0.0
        self._int_theta = 0.0
        self._int_dtheta = 0.0

        # Simple derivative filtering for alpha (on measurement)
        self._dalpha_filt = 0.0
        self._dalpha_tau = 0.02  # [s]

    def reset(self) -> None:
        self._int_alpha = 0.0
        self._int_theta = 0.0
        self._int_dtheta = 0.0
        self._dalpha_filt = 0.0

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

        This uses:
          - Arm: cascade PI (theta) -> dtheta_ref, then PI (dtheta)
          - Pendulum: PID (typically PD is enough; ki_alpha defaults to 0)
        """

        Ts = self.Ts
        if Ts <= 0.0:
            Ts = 1e-3

        # --- Errors ---
        e_theta = _wrap_to_pi(self.theta_ref_rad - float(theta))
        e_alpha = _wrap_to_pi(self.alpha_ref_rad - float(alpha))

        # --- Filtered derivative for pendulum (use measurement derivative) ---
        # Low-pass filter: y_dot_filt = (tau/(tau+Ts))*y_dot_filt + (Ts/(tau+Ts))*y_dot
        tau = self._dalpha_tau
        if tau > 0.0:
            a = tau / (tau + Ts)
            b = Ts / (tau + Ts)
            self._dalpha_filt = a * self._dalpha_filt + b * float(dalpha)
            dalpha_used = self._dalpha_filt
        else:
            dalpha_used = float(dalpha)

        # --- Outer loop: theta -> desired dtheta ---
        # PI (no derivative) to avoid double-derivative when inner loop uses velocity.
        dtheta_ref = self.kp_theta_outer * e_theta + \
            self.ki_theta_outer * self._int_theta

        # --- Inner loop: dtheta -> voltage contribution ---
        e_dtheta = float(dtheta_ref) - float(dtheta)
        v_theta = self.kp_dtheta_inner * e_dtheta + \
            self.ki_dtheta_inner * self._int_dtheta

        # --- Pendulum stabilization voltage contribution ---
        # Derivative term uses -dalpha (derivative of error when alpha_ref is constant)
        v_alpha = (
            self.kp_alpha * e_alpha
            + self.ki_alpha * self._int_alpha
            - self.kd_alpha * dalpha_used
        )

        v_unsat = float(v_theta + v_alpha)

        # --- Saturation ---
        v_sat = max(-self.v_limit, min(self.v_limit, v_unsat))

        # --- Anti-windup (conditional integration) ---
        # If saturated and the integrator would push further into saturation, freeze that integrator.
        def _should_integrate(v_cmd: float, v_cmd_sat: float, err: float) -> bool:
            if v_cmd == v_cmd_sat:
                return True
            # saturated high and error wants more positive -> don't integrate
            if v_cmd_sat >= self.v_limit - 1e-12 and err > 0.0:
                return False
            # saturated low and error wants more negative -> don't integrate
            if v_cmd_sat <= -self.v_limit + 1e-12 and err < 0.0:
                return False
            return True

        if _should_integrate(v_unsat, v_sat, e_alpha):
            self._int_alpha += e_alpha * Ts
        if _should_integrate(v_unsat, v_sat, e_theta):
            self._int_theta += e_theta * Ts
        if _should_integrate(v_unsat, v_sat, e_dtheta):
            self._int_dtheta += e_dtheta * Ts

        return v_sat
