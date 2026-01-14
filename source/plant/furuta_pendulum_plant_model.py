from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[2]))

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import dill

NPZ_FILE_PATH = Path(__file__).parent / "furuta_pendulum_plant_model.npz"


class FurutaPendulum:
    """
    Furuta pendulum nonlinear dynamics derived by Lagrange (symbolic),
    including DC motor electrical dynamics (R-L + back-EMF) and torque generation.

    State x = [theta, alpha, theta_dot, alpha_dot, i]
    Input u = v (motor voltage)
    """

    def __init__(self, params: dict):
        self.p = dict(params)

        self.input_time_series = None
        self.input_value_series = None

        self.f_thdd = None
        self.f_aldd = None
        self.f_di = None

        self._param_tuple = None

        if not NPZ_FILE_PATH.exists():
            self._build_symbolic_model()

            pickled = dill.dumps(self)
            np.savez_compressed(NPZ_FILE_PATH,
                                pickled_instance=np.array([pickled], dtype=object))

        else:
            data = np.load(NPZ_FILE_PATH, allow_pickle=True)
            loaded_instance: FurutaPendulum = dill.loads(
                data["pickled_instance"][0])
            self.f_thdd = loaded_instance.f_thdd
            self.f_aldd = loaded_instance.f_aldd
            self.f_di = loaded_instance.f_di
            self._param_tuple = loaded_instance._param_tuple

    def _build_symbolic_model(self):
        p = self.p

        # Symbols
        theta, alpha = sp.symbols("theta alpha", real=True)
        thd, ald = sp.symbols("thd ald", real=True)
        thdd, aldd = sp.symbols("thdd aldd", real=True)
        i, v = sp.symbols("i v", real=True)

        # Parameters (symbols)
        m_p, L_r, L_p, g = sp.symbols(
            "m_p L_r L_p g", positive=True, real=True)
        J_r, J_m, n = sp.symbols("J_r J_m n", positive=True, real=True)
        rod_rad, p_rad = sp.symbols("rod_rad p_rad", positive=True, real=True)

        D_r, D_p, mu_m = sp.symbols(
            "D_r D_p mu_m", positive=True, real=True)
        R_m, L_m, K_t, K_b = sp.symbols(
            "R_m L_m K_t K_b", positive=True, real=True)

        # Useful
        l_p = L_p / 2

        # Unit vectors in inertial frame
        e_r = sp.Matrix([sp.cos(theta), sp.sin(theta), 0]
                        )              # arm direction
        e_th = sp.Matrix([-sp.sin(theta), sp.cos(theta), 0]
                         )            # azimuth tangential
        e_z = sp.Matrix([0, 0, 1])

        # Pendulum rod direction (alpha=0 is upright +z)
        # Use rotation about the arm direction `e_r` so the pendulum
        # pitches around the arm axis (matches Furuta convention).
        # Rotation of +alpha moves rod toward -e_th (see derivation):
        # u_hat = cos(alpha)*e_z - sin(alpha)*e_th
        u_hat = sp.cos(alpha) * e_z - sp.sin(alpha) * e_th

        # Pivot at arm end, COM position
        r_piv = L_r * e_r
        r_c = r_piv + l_p * u_hat

        # Velocity of COM
        v_c = sp.diff(r_c, theta) * thd + sp.diff(r_c, alpha) * ald

        # Angular velocity of pendulum body:
        # arm yaw contributes thd about z, pendulum pitch contributes ald about e_r
        omega = thd * e_z + ald * e_r

        # Inertia about pendulum COM:
        # I_perp_cm matches your formula part: m*(r^2/4 + L^2/12)
        I_perp_cm = m_p * (p_rad**2 / 4 + L_p**2 / 12)
        # Approx rod-axis inertia (solid cylinder): 1/2 m r^2 (small but included)
        I_long = m_p * (p_rad**2 / 2)

        omega_sq = (omega.dot(omega))
        omega_par = (omega.dot(u_hat))
        omega_perp_sq = sp.simplify(omega_sq - omega_par**2)

        T_p_trans = sp.Rational(1, 2) * m_p * (v_c.dot(v_c))
        T_p_rot = sp.Rational(1, 2) * I_perp_cm * omega_perp_sq + \
            sp.Rational(1, 2) * I_long * omega_par**2

        # Arm kinetic energy about base (vertical axis): include arm rod + reflected motor inertia
        Jm_eq = n**2 * J_m
        T_arm = sp.Rational(1, 2) * (J_r + Jm_eq) * thd**2

        T = sp.simplify(T_arm + T_p_trans + T_p_rot)

        # Potential energy (z is upward)
        V = sp.simplify(m_p * g * r_c[2])

        Lagr = sp.simplify(T - V)

        # Rayleigh dissipation (viscous)
        Dth = D_r + n**2 * mu_m
        R = sp.Rational(1, 2) * Dth * thd**2 + sp.Rational(1, 2) * D_p * ald**2

        # Generalized forces (input torque on theta only)
        tau = n * K_t * i
        Q = sp.Matrix([tau, 0])

        q = sp.Matrix([theta, alpha])
        qd = sp.Matrix([thd, ald])
        qdd = sp.Matrix([thdd, aldd])

        # Lagrange equations with Rayleigh:
        # d/dt(dL/dqd) - dL/dq + dR/dqd = Q
        dL_dq = sp.Matrix([sp.diff(Lagr, theta), sp.diff(Lagr, alpha)])
        dL_dqd = sp.Matrix([sp.diff(Lagr, thd), sp.diff(Lagr, ald)])
        dR_dqd = sp.Matrix([sp.diff(R, thd), sp.diff(R, ald)])

        # time derivative of dL/dqd using chain rule
        ddt_dL_dqd = sp.Matrix([
            sp.diff(dL_dqd[k], theta) * thd
            + sp.diff(dL_dqd[k], alpha) * ald
            + sp.diff(dL_dqd[k], thd) * thdd
            + sp.diff(dL_dqd[k], ald) * aldd
            for k in range(2)
        ])

        eom = sp.simplify(ddt_dL_dqd - dL_dq + dR_dqd - Q)  # = 0

        # Extract M(q) and rhs so that M*qdd = rhs
        M = sp.simplify(eom.jacobian(qdd))
        h = sp.simplify(eom - M * qdd)
        rhs = sp.simplify(-h)

        # Solve for qdd
        qdd_sol = sp.simplify(M.LUsolve(rhs))
        thdd_sol = sp.simplify(qdd_sol[0])
        aldd_sol = sp.simplify(qdd_sol[1])

        # Motor electrical dynamics: v = R i + L di/dt + Kb * omega_m ; omega_m = n*thd
        di_sol = sp.simplify((v - R_m * i - K_b * (n * thd)) / L_m)

        # Lambdify
        sym_list = [
            theta, alpha, thd, ald, i, v,
            # params
            m_p, L_r, L_p, g, J_r, J_m, n, rod_rad, p_rad,
            D_r, D_p, mu_m, R_m, L_m, K_t, K_b
        ]

        self.f_thdd = sp.lambdify(sym_list, thdd_sol, "numpy")
        self.f_aldd = sp.lambdify(sym_list, aldd_sol, "numpy")
        self.f_di = sp.lambdify(sym_list, di_sol, "numpy")

        # Cache numeric parameter tuple order
        self._param_tuple = (
            float(p["m_p"]), float(p["L_r"]), float(p["L_p"]), float(p["g"]),
            float(p["J_r"]), float(p["J_m"]), float(
                p["n"]), float(p["rod_rad"]), float(p["p_rad"]),
            float(p["D_r"]), float(p["D_p"]), float(
                p["mu_m"]), float(p["R_m"]),
            float(p["L_m"]), float(p["K_t"]), float(p["K_b"])
        )

    def dynamics(self, t, x, v_in):
        """
        Continuous-time dynamics.
        x = [theta, alpha, theta_dot, alpha_dot, i]
        """
        theta, alpha, thd, ald, i = x
        v = float(v_in)

        (m_p, L_r, L_p, g,
         J_r, J_m, n, rod_rad, p_rad,
         D_r, D_p, mu_m, R_m,
         L_m, K_t, K_b) = self._param_tuple

        thdd = self.f_thdd(theta, alpha, thd, ald, i, v,
                           m_p, L_r, L_p, g, J_r, J_m, n, rod_rad, p_rad,
                           D_r, D_p, mu_m, R_m, L_m, K_t, K_b)

        aldd = self.f_aldd(theta, alpha, thd, ald, i, v,
                           m_p, L_r, L_p, g, J_r, J_m, n, rod_rad, p_rad,
                           D_r, D_p, mu_m, R_m, L_m, K_t, K_b)

        di = self.f_di(theta, alpha, thd, ald, i, v,
                       m_p, L_r, L_p, g, J_r, J_m, n, rod_rad, p_rad,
                       D_r, D_p, mu_m, R_m, L_m, K_t, K_b)

        return np.array([thd, ald, float(thdd), float(aldd), float(di)], dtype=float)

    def simulate(
        self,
        x0,
        t_span,
        v_func,
        v_func_time_step,
        dt=0.001,
        rtol=1e-7,
        atol=1e-9
    ):
        """
        Simulate with input voltage v_func(t, x)->v.
        Returns t, X arrays.
        """
        v_func_time_count = 0.0
        input_time_series = []
        input_value_series = []
        v_saturated = 0.0

        def f(t, x):
            nonlocal v_func_time_count
            nonlocal input_time_series
            nonlocal input_value_series
            nonlocal v_saturated

            if t >= v_func_time_count:

                v = v_func(t, x)
                v_saturated = max(self.p["V_min"], min(self.p["V_max"], v))

                input_time_series.append(t)
                input_value_series.append(v_saturated)

                v_func_time_count += v_func_time_step

            return self.dynamics(t, x, v_saturated)

        sol = solve_ivp(
            f, t_span, np.array(x0, dtype=float),
            method="RK45", max_step=dt, rtol=rtol, atol=atol
        )

        self.input_time_series = np.array(input_time_series, dtype=float)
        self.input_value_series = np.array(input_value_series, dtype=float)

        return sol.t, sol.y.T


# =========================
# シミュレーションパラメータ
# =========================
params = dict(
    R_m=21.7,
    mu_m=3.08e-6,
    K_t=0.042,
    K_b=0.182,
    J_m=4e-6,
    L_m=4.98e-3,

    rod_rad=0.003,
    p_rad=0.0045,
    L_p=0.126,
    L_r=0.103,
    m_p=0.024,
    m_r=0.095,

    J_r=0.095 * (0.003**2 / 4 + 0.103**2 / 12) + 0.095 * (0.103 / 2)**2,

    D_r=1.88e-4,
    D_p=8e-6,

    n=1,
    g=9.81,

    V_max=12.0,
    V_min=-12.0,
)

# Sampled-data controller helper: executes the provided control law only at
# discrete sample instants and holds the value (zero-order hold) between
# samples. This ensures the controller is evaluated at a fixed period (Ts)
# while solve_ivp may request continuous-time dynamics at arbitrary times.


class SampledController:
    def __init__(self, Ts, control_law, sat=None):
        self.Ts = float(Ts)
        self.control_law = control_law
        self.last_idx = None
        self.last_u = 0.0
        self.sat = sat

    def __call__(self, t, x):
        # determine which sample interval this time belongs to
        idx = int(np.floor(t / self.Ts + 1e-12))
        # if a new sample instant, evaluate control law with current state
        if self.last_idx is None or idx != self.last_idx:
            # ensure x is a numpy array for indexing
            x_arr = np.array(x, dtype=float)
            self.last_u = float(self.control_law(t, x_arr))
            # optional saturation
            if self.sat is not None:
                lo, hi = self.sat
                self.last_u = max(lo, min(hi, self.last_u))
            self.last_idx = idx
        # hold the last computed value until next sample
        return self.last_u
