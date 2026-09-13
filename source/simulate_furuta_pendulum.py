from __future__ import annotations

import os
import sys
import threading
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from external_libraries.simulation_manager.visualize.simulation_plotter_dash import SimulationPlotterDash
from source.plot_pendulum_move import plot_matplotlib, plot_plotly

from source.plant.furuta_pendulum_plant_model import (
    FurutaPendulum,
    params,
    SampledController,
)

from source.controller.furuta_pendulum_pid_controller import FurutaPendulum_PID_Controller

SIMULATION_TIME_STEP = 0.001  # シミュレーションの時間刻み幅（秒）
SIMULATION_END_TIME = 10.0    # シミュレーションの終了時間（秒）
PLAYBACK_FPS = 1000  # 3Dプロット再生時のフレームレート（FPS）

CONTROLLER_TIME_STEP = 0.005  # コントローラの時間刻み幅（秒）

FAST_RESTART = False  # 物理モデルの再構築をせず、前回の結果を使用する場合はTrueに設定
SIL_MODE = False  # C++コードをSIL検証する場合はTrueに設定
SIL_TOLERANCE = 1e-5  # SIL検証の許容誤差

# 物理モデル
print("Building symbolic physics model...")
model = FurutaPendulum(params, fast_restart=FAST_RESTART)
print("Model built successfully.")

# 初期値
# x = [theta, alpha, theta_dot, alpha_dot, i]
x0 = [0.0, np.deg2rad(10.0), 0.0, 0.0, 0.0]

# コントローラー
controller = FurutaPendulum_PID_Controller(Ts=CONTROLLER_TIME_STEP)

if SIL_MODE:
    from external_libraries.MCAP_Playground.helper.SIL.SIL_operator import SIL_Operator

    current_dir = os.path.dirname(__file__)
    generator = SIL_Operator("furuta_pendulum_pid_controller.py", current_dir)
    generator.build_SIL_code(build_type="Debug")

    import FurutaPendulumPidControllerSIL
    FurutaPendulumPidControllerSIL.initialize(CONTROLLER_TIME_STEP)

    voltage_cpp = 0.0

# プロット
plotter = SimulationPlotterDash()


def feedback_law(t, x):

    theta = x[0]
    alpha = x[1]
    dtheta = x[2]
    dalpha = x[3]

    theta_ref = 0.0

    if t >= 5.0:
        theta_ref = np.deg2rad(45.0)

    controller.set_theta_reference_rad(theta_ref)
    if SIL_MODE:
        FurutaPendulumPidControllerSIL.set_theta_reference_rad(theta_ref)

    voltage = controller.calculate_manipulation(
        theta, alpha, dtheta, dalpha)
    if SIL_MODE:
        voltage_cpp = FurutaPendulumPidControllerSIL.calculate_manipulation(
            theta, alpha, dtheta, dalpha)

        assert voltage == pytest.approx(
            voltage_cpp, abs=SIL_TOLERANCE), "Voltage mismatch"

        plotter.append_name(voltage_cpp, "voltage_cpp")

    return voltage


sampled_controller = SampledController(
    CONTROLLER_TIME_STEP, feedback_law, sat=(-12.0, 12.0))

print("Running simulation...")
t_sim, X_sim = model.simulate(
    x0=x0,
    t_span=(0.0, SIMULATION_END_TIME),
    v_func=sampled_controller,
    v_func_time_step=sampled_controller.Ts,
    dt=SIMULATION_TIME_STEP)
print(f"Simulation complete. {len(t_sim)} time steps.")

# シミュレーション結果から theta, alpha を抽出
time_series = t_sim
theta = X_sim[:, 0]  # arm angle
alpha = X_sim[:, 1]  # pendulum angle
voltage = model.input_value_series
voltage_time = model.input_time_series

# Quick check prints
print("Final state:", X_sim[-1])
print("Max |alpha| [deg]:", np.rad2deg(np.max(np.abs(alpha))))

# 可視化用のアーム・振子長さ（物理パラメータから取得）
L_arm = params["L_r"]
L_pend = params["L_p"]

# 波形表示


plotter.append_sequence_name(theta, "theta")
plotter.append_sequence_name(alpha, "alpha")
plotter.append_sequence_name(voltage, "voltage")

plotter.assign("theta", row=0, column=0, position=(0, 0),
               x_sequence=time_series, label="theta")
plotter.assign("alpha", row=0, column=0, position=(1, 0),
               x_sequence=time_series, label="alpha")
plotter.assign("voltage", row=0, column=0, position=(2, 0),
               x_sequence=voltage_time, label="voltage")
if SIL_MODE:
    plotter.assign("voltage_cpp", row=0, column=0, position=(2, 0),
                   x_sequence=voltage_time, label="voltage_cpp")

plotter.plot(suptitle="Furuta Pendulum Simulation Results")

plot_plotly(time_series, theta, alpha, L_arm, L_pend)
# plot_matplotlib(time_series, theta, alpha, L_arm, L_pend)
