"""
Furuta Pendulum 3D Visualization Module

Provides two visualization backends:
  - matplotlib: interactive 3D animation with Play/Stop buttons and a time slider
  - plotly:     standalone HTML with a time slider (drag to rotate camera)
"""
from __future__ import annotations

import json
import os
import webbrowser
from datetime import datetime
from pathlib import Path

import numpy as np

PLAYBACK_FPS = 1000  # default playback frame rate for both backends

# ---------------------------------------------------------------------------
# Common geometry helper
# ---------------------------------------------------------------------------


def _compute_geometry(theta: np.ndarray, alpha: np.ndarray,
                      L_arm: float, L_pend: float):
    """Return (arm_tip_xyz, pend_tip_xyz) for every time step."""
    x_arm = L_arm * np.cos(theta)
    y_arm = L_arm * np.sin(theta)
    z_arm = np.zeros_like(theta)

    x_p = x_arm + L_pend * np.sin(alpha) * np.sin(theta)
    y_p = y_arm - L_pend * np.sin(alpha) * np.cos(theta)
    z_p = L_pend * np.cos(alpha)

    return (x_arm, y_arm, z_arm), (x_p, y_p, z_p)


# ===========================================================================
# Backend 1: matplotlib
# ===========================================================================

def plot_matplotlib(time_series: np.ndarray,
                    theta: np.ndarray,
                    alpha: np.ndarray,
                    L_arm: float,
                    L_pend: float,
                    playback_fps: int = PLAYBACK_FPS) -> None:
    """
    Show a matplotlib 3D animation of the Furuta pendulum.

    Controls:
        Play  – start playback
        Stop  – pause playback
        Slider – seek to any time step
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button, Slider

    (x_arm, y_arm, z_arm), (x_p, y_p, z_p) = _compute_geometry(
        theta, alpha, L_arm, L_pend)

    N = len(time_series)
    dt_sample = (time_series[1] - time_series[0]) if N > 1 else 1e-3

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    plt.subplots_adjust(bottom=0.32)

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
    time_text = ax.text2D(
        0.02, 0.95, "",
        transform=ax.transAxes, fontsize=12,
        bbox=dict(facecolor="white", alpha=0.7),
    )

    idx = 0
    running = False

    def draw(i: int):
        arm_line.set_data([0, x_arm[i]], [0, y_arm[i]])
        arm_line.set_3d_properties([0, z_arm[i]])
        pend_line.set_data([x_arm[i], x_p[i]], [y_arm[i], y_p[i]])
        pend_line.set_3d_properties([z_arm[i], z_p[i]])
        time_text.set_text(f"t = {time_series[i]:.2f} [s]")
        fig.canvas.draw_idle()

    ax_play = plt.axes([0.35, 0.14, 0.12, 0.06])
    ax_stop = plt.axes([0.53, 0.14, 0.12, 0.06])
    btn_play = Button(ax_play, "Play")
    btn_stop = Button(ax_stop, "Stop")

    def play(_):
        nonlocal running
        running = True

    def stop(_):
        nonlocal running
        running = False

    btn_play.on_clicked(play)
    btn_stop.on_clicked(stop)

    ax_slider_time = plt.axes([0.20, 0.20, 0.60, 0.03])
    slider_time = Slider(ax_slider_time, "Step", 0, N - 1,
                         valinit=0, valstep=1)

    def on_time_slider(val):
        nonlocal idx
        idx = int(val)
        draw(idx)

    slider_time.on_changed(on_time_slider)
    draw(0)

    playback_step = max(1, int(round((1.0 / playback_fps) / dt_sample)))

    def update_timer():
        nonlocal idx
        if running:
            idx = min(idx + playback_step, N - 1)
            slider_time.set_val(idx)

    timer = fig.canvas.new_timer(interval=int(1000.0 / playback_fps))
    timer.add_callback(update_timer)
    timer.start()

    plt.show()


# ===========================================================================
# Backend 2: Plotly (standalone HTML)
# ===========================================================================

_CACHE_FOLDER = Path.cwd() / "cache"

_PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.27.0.min.js"


def plot_plotly(time_series: np.ndarray,
                theta: np.ndarray,
                alpha: np.ndarray,
                L_arm: float,
                L_pend: float,
                title: str = "Furuta Pendulum 3D",
                output_dir: Path | None = None) -> None:
    """
    Generate a standalone HTML file with an interactive Plotly 3D view.

    Features:
      - Drag to rotate camera (built-in Plotly 3D behaviour)
      - Slider to seek any time step
      - XYZ axis grid lines identical to the matplotlib view
    """
    (x_arm_arr, y_arm_arr, z_arm_arr), (x_p_arr, y_p_arr, z_p_arr) = \
        _compute_geometry(theta, alpha, L_arm, L_pend)

    N = len(time_series)
    axis_limit = (L_arm + L_pend) / 2.0 * 1.5

    # Subsample the stored frames to keep the HTML size reasonable.
    # Keep at most 2000 frames; for shorter sims keep every frame.
    max_frames = 2000
    step = max(1, N // max_frames)
    indices = list(range(0, N, step))
    if indices[-1] != N - 1:
        indices.append(N - 1)

    # Build coordinate arrays for each sampled frame.
    # Each frame stores: arm [ox,tip], pend [armtip, pendtip]
    frames_data = []
    for i in indices:
        frames_data.append({
            "t": float(time_series[i]),
            "ax": [0.0, float(x_arm_arr[i])],
            "ay": [0.0, float(y_arm_arr[i])],
            "az": [0.0, float(z_arm_arr[i])],
            "px": [float(x_arm_arr[i]), float(x_p_arr[i])],
            "py": [float(y_arm_arr[i]), float(y_p_arr[i])],
            "pz": [float(z_arm_arr[i]), float(z_p_arr[i])],
        })

    frames_json = json.dumps(frames_data)

    axis_cfg = dict(
        range=[-axis_limit, axis_limit],
        showgrid=True,
        gridcolor="#cccccc",
        zeroline=True,
        zerolinecolor="#888888",
        showline=True,
        linecolor="#333333",
        showspikes=False,
        backgroundcolor="white",
        showbackground=True,
    )

    scene_layout = dict(
        xaxis=dict(**axis_cfg, title="X"),
        yaxis=dict(**axis_cfg, title="Y"),
        zaxis=dict(**axis_cfg, title="Z"),
        bgcolor="white",
        camera=dict(eye=dict(x=1.4, y=1.4, z=0.8)),
        aspectmode="cube",
    )

    initial = frames_data[0]

    import plotly.graph_objects as go

    fig = go.Figure(data=[
        go.Scatter3d(
            x=initial["ax"], y=initial["ay"], z=initial["az"],
            mode="lines",
            line=dict(color="red", width=8),
            name="Arm",
        ),
        go.Scatter3d(
            x=initial["px"], y=initial["py"], z=initial["pz"],
            mode="lines",
            line=dict(color="blue", width=8),
            name="Pendulum",
        ),
        # Invisible origin marker so the scene range is anchored
        go.Scatter3d(
            x=[0], y=[0], z=[0],
            mode="markers",
            marker=dict(size=3, color="black"),
            name="Origin",
            showlegend=False,
        ),
    ])

    fig.update_layout(
        title=title,
        scene=scene_layout,
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=650,
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.8)",
                    bordercolor="black", borderwidth=1),
        margin=dict(l=0, r=0, t=40, b=0),
    )

    fig_json = fig.to_json()

    # -----------------------------------------------------------------------
    # Standalone HTML with embedded JS slider
    # -----------------------------------------------------------------------
    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{title}</title>
  <script src="{_PLOTLY_CDN}"></script>
  <style>
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: sans-serif; background: white; }}
    #controls {{
      display: flex;
      align-items: center;
      gap: 16px;
      padding: 8px 16px;
      background: #f5f5f5;
      border-bottom: 1px solid #ddd;
    }}
    #controls label {{ font-size: 14px; }}
    #controls button {{
      padding: 4px 14px;
      font-size: 14px;
      cursor: pointer;
    }}
    #time-slider {{ flex: 1; }}
    #time-display {{ font-size: 13px; font-family: monospace; min-width: 90px; }}
  </style>
</head>
<body>
  <div id="controls">
    <button id="btn-play">Play</button>
    <button id="btn-stop">Stop</button>
    <label for="time-slider">Time:</label>
    <input id="time-slider" type="range" min="0" max="{len(frames_data) - 1}"
           value="0" step="1">
    <span id="time-display">t = 0.000 s</span>
  </div>
  <div id="graph" style="height:650px;"></div>

  <script>
  (function() {{
    var framesData = {frames_json};
    var figData = {fig_json};

    var gd = document.getElementById("graph");
    Plotly.newPlot(gd, figData.data, figData.layout,
      {{scrollZoom: true, displayModeBar: true}});

    var slider = document.getElementById("time-slider");
    var display = document.getElementById("time-display");
    var btnPlay = document.getElementById("btn-play");
    var btnStop = document.getElementById("btn-stop");

    function updateFrame(idx) {{
      var f = framesData[idx];
      display.textContent = "t = " + f.t.toFixed(3) + " s";
      slider.value = idx;
      Plotly.restyle(gd, {{
        x: [f.ax, f.px, [0]],
        y: [f.ay, f.py, [0]],
        z: [f.az, f.pz, [0]],
      }}, [0, 1, 2]);
    }}

    slider.addEventListener("input", function() {{
      updateFrame(parseInt(slider.value, 10));
    }});

    // Real-time playback: advance frames based on wall-clock elapsed time
    var rafId = null;
    var playStartWall = 0;
    var playStartSimTime = 0;
    var currentIdx = 0;

    function findFrameIdx(simTime) {{
      // Binary search for the frame closest to simTime
      var lo = 0, hi = framesData.length - 1;
      while (lo < hi) {{
        var mid = (lo + hi + 1) >> 1;
        if (framesData[mid].t <= simTime) lo = mid; else hi = mid - 1;
      }}
      return lo;
    }}

    function playLoop(wallNow) {{
      var elapsed = (wallNow - playStartWall) / 1000.0; // seconds
      var targetSimTime = playStartSimTime + elapsed;
      var idx = findFrameIdx(targetSimTime);
      if (idx !== currentIdx) {{
        currentIdx = idx;
        updateFrame(idx);
      }}
      if (idx < framesData.length - 1) {{
        rafId = requestAnimationFrame(playLoop);
      }} else {{
        rafId = null;
      }}
    }}

    btnPlay.addEventListener("click", function() {{
      if (rafId !== null) return;
      currentIdx = parseInt(slider.value, 10);
      // If at the end, restart from beginning
      if (currentIdx >= framesData.length - 1) currentIdx = 0;
      playStartWall = performance.now();
      playStartSimTime = framesData[currentIdx].t;
      rafId = requestAnimationFrame(playLoop);
    }});

    btnStop.addEventListener("click", function() {{
      if (rafId !== null) {{
        cancelAnimationFrame(rafId);
        rafId = null;
      }}
    }});
  }})();
  </script>
</body>
</html>
"""

    if output_dir is None:
        output_dir = _CACHE_FOLDER
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = output_dir / f"{timestamp}_{title.replace(' ', '_')}.html"
    out_path.write_text(html, encoding="utf-8")
    print(f"Plotly 3D HTML saved to: {out_path}")

    if not os.environ.get("CI") and not os.environ.get("GITHUB_ACTIONS"):
        if not webbrowser.open(out_path.as_uri()):
            print(
                f"Warning: Could not open browser. Open manually: {out_path}")
