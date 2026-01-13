from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

from external_libraries.simulation_manager.visualize.simulation_plotter import SimulationPlotter

from source.plant.furuta_pendulum_plant_model import (
    FurutaPendulum,
    params,
    SampledController,
)

from source.controller.furuta_pendulum_pid_controller import FurutaPendulum_PID_Controller

SIMULATION_TIME_STEP = 0.001  # シミュレーションの時間刻み幅（秒）
SIMULATION_END_TIME = 5.0    # シミュレーションの終了時間（秒）
PLAYBACK_FPS = 200  # 3Dプロット再生時のフレームレート（FPS）

CONTROLLER_TIME_STEP = 0.005  # コントローラの時間刻み幅（秒）

# =========================
# 物理モデルの構築とシミュレーション実行
# =========================
print("Building symbolic physics model...")
model = FurutaPendulum(params, use_alt_backemf=True,
                       use_alt_arm_damping=True)
print("Model built successfully.")

# Initial condition:
# x = [theta, alpha, theta_dot, alpha_dot, i]
x0 = [0.0, np.deg2rad(10.0), 0.0, 0.0, 0.0]

# コントローラー
controller = FurutaPendulum_PID_Controller()


def feedback_law(t, x):

    theta = x[0]
    alpha = x[1]
    dtheta = x[2]
    dalpha = x[3]

    voltage = controller.calculate_manipulation(
        theta, alpha, dtheta, dalpha)

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

N = len(time_series)
dt = time_series[1] - \
    time_series[0] if len(time_series) > 1 else SIMULATION_TIME_STEP
# Quick check prints
print("Final state:", X_sim[-1])
print("Max |alpha| [deg]:", np.rad2deg(np.max(np.abs(alpha))))

# 可視化用のアーム・振子長さ（物理パラメータから取得）
L_arm = params["L_r"]
L_pend = params["L_p"]

# =========================
# Figure & Axes
# =========================
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection="3d")

plt.subplots_adjust(bottom=0.32)

# 軸の範囲を物理パラメータに基づいて設定
axis_limit = (L_arm + L_pend) / 2.0 * 1.5
ax.set_xlim([-axis_limit, axis_limit])
ax.set_ylim([-axis_limit, axis_limit])
ax.set_zlim([-axis_limit, axis_limit])
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.view_init(elev=25, azim=45)

arm_line, = ax.plot([], [], [], "r-", lw=3)
pend_line, = ax.plot([], [], [], "b-", lw=3)

# =========================
# 時刻表示（画面固定）
# =========================
time_text = ax.text2D(
    0.02, 0.95,
    "",
    transform=ax.transAxes,
    fontsize=12,
    bbox=dict(facecolor="white", alpha=0.7)
)

# =========================
# 状態管理
# =========================
idx = 0
running = False


# =========================
# 描画更新
# =========================
def draw(i):
    th = theta[i]
    al = alpha[i]

    x_arm = L_arm * np.cos(th)
    y_arm = L_arm * np.sin(th)
    z_arm = 0.0

    # Pendulum tip computed using the same kinematic convention used
    # in the symbolic model: u_hat = cos(alpha)*e_z - sin(alpha)*e_th
    # which expands to components:
    # [ sin(alpha)*sin(th), -sin(alpha)*cos(th), cos(alpha) ]
    x_p = x_arm + L_pend * (np.sin(al) * np.sin(th))
    y_p = y_arm - L_pend * (np.sin(al) * np.cos(th))
    z_p = L_pend * np.cos(al)

    arm_line.set_data([0, x_arm], [0, y_arm])
    arm_line.set_3d_properties([0, z_arm])

    pend_line.set_data([x_arm, x_p], [y_arm, y_p])
    pend_line.set_3d_properties([z_arm, z_p])

    time_text.set_text(f"t = {time_series[i]:.2f} [s]")

    fig.canvas.draw_idle()

# =========================
# タイマー処理（再生はPLAYBACK_FPSで行う）
# =========================


def update_timer():
    global idx
    if running:
        idx = min(idx + playback_step, N - 1)
        slider_time.set_val(idx)

# タイマーは後で作成（slider_timeが定義された後）


# =========================
# 再生 / 停止
# =========================
# Buttons moved below the slider (centered under slider area)
ax_play = plt.axes([0.35, 0.14, 0.12, 0.06])
ax_stop = plt.axes([0.53, 0.14, 0.12, 0.06])

btn_play = Button(ax_play, "Play")
btn_stop = Button(ax_stop, "Stop")


def play(event):
    global running
    running = True


def stop(event):
    global running
    running = False


btn_play.on_clicked(play)
btn_stop.on_clicked(stop)

# =========================
# 時間シークバー
# =========================
# Slider placed above the buttons and stretched wider
ax_slider_time = plt.axes([0.20, 0.20, 0.60, 0.03])
slider_time = Slider(
    ax_slider_time,
    "Step",
    0,
    N - 1,
    valinit=0,
    valstep=1
)


def on_time_slider(val):
    global idx
    idx = int(val)
    draw(idx)


slider_time.on_changed(on_time_slider)

# 3D描画
draw(0)
dt_sample = dt if N > 1 else SIMULATION_TIME_STEP
playback_step = max(1, int(round((1.0 / PLAYBACK_FPS) / dt_sample)))

timer = fig.canvas.new_timer(interval=int(1000.0 / PLAYBACK_FPS))
timer.add_callback(update_timer)
timer.start()

# 波形表示
plotter = SimulationPlotter()

plotter.append_sequence_name(theta, "theta")
plotter.append_sequence_name(alpha, "alpha")
plotter.append_sequence_name(voltage, "voltage")

plotter.assign("theta", column=0, row=0, position=(0, 0),
               x_sequence=time_series, label="theta")
plotter.assign("alpha", column=0, row=0, position=(1, 0),
               x_sequence=time_series, label="alpha")
plotter.assign("voltage", column=0, row=0, position=(2, 0),
               x_sequence=voltage_time, label="voltage")

plotter.pre_plot(suptitle="Furuta Pendulum Simulation Results")

plt.show()
