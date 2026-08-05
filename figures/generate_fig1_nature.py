"""
Figure 1 (PDF) – Simulation setup for spectral DCM / calcium imaging paper.

2-row layout:
  Row 1: A  4-node network diagram (ground-truth network)
         B  PSD overlay — neural state (solid) + Ca²⁺ (dashed, darker)
  Row 2: C  Time series overlay — neural state (lighter) + Ca²⁺ (darker, same offset)

Save: figures/Fig1.png
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D

sys.path.insert(0, "/Users/hjpark/Dropbox/matlabwork/mnet0.92/dcmcai")
import nature_style as ns

OUT_DIR = "/Users/hjpark/Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures"
os.makedirs(OUT_DIR, exist_ok=True)

# ── True connectivity matrix ──────────────────────────────────────────────────
Atrue = np.array([[-0.5,  0.0, -0.3, -0.1],
                  [ 0.4, -0.5,  0.2,  0.0],
                  [ 0.0,  0.2, -0.5, -0.1],
                  [ 0.1,  0.3,  0.0, -0.5]])

# ── Simulation ────────────────────────────────────────────────────────────────
np.random.seed(42)
TR     = 0.5
T      = 2500
sigma  = 0.05
tau_ca = 0.75
N      = 4

x = np.zeros((T, N))
for t in range(1, T):
    dW = np.random.randn(N)
    x[t] = x[t-1] + (Atrue @ x[t-1]) * TR + sigma * np.sqrt(TR) * dW

h_decay = np.exp(-TR / tau_ca)
y = np.zeros_like(x)
for t in range(1, T):
    y[t] = h_decay * y[t-1] + (1 - h_decay) * x[t]

x_std = x.std(axis=0)
x_z = (x - x.mean(axis=0)) / x_std
y_z = (y - y.mean(axis=0)) / x_std

# ── Analytical PSD ────────────────────────────────────────────────────────────
def _analytical_psd(freqs_hz):
    Q = sigma**2 * np.eye(N)
    psd_x = np.zeros((N, len(freqs_hz)))
    psd_y = np.zeros((N, len(freqs_hz)))
    for fi, f in enumerate(freqs_hz):
        w = 2.0 * np.pi * f
        G = np.linalg.solve(1j * w * np.eye(N) - Atrue, np.eye(N))
        Sxx = G @ Q @ G.conj().T
        H2 = 1.0 / (1.0 + (w * tau_ca) ** 2)
        for i in range(N):
            psd_x[i, fi] = Sxx[i, i].real
            psd_y[i, fi] = H2 * Sxx[i, i].real
    return psd_x, psd_y

_freqs = np.linspace(0.0078125, 0.25, 64)
_psd_x, _psd_y = _analytical_psd(_freqs)

# ── Colors ────────────────────────────────────────────────────────────────────
ns.apply()

node_labels = ['R1', 'R2', 'R3', 'R4']
colors4 = [ns.COLORS['blue'], ns.COLORS['orange'],
           ns.COLORS['green'], ns.COLORS['vermil']]

def darken(color, factor=0.55):
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)

dark4 = [darken(c) for c in colors4]

node_pos = {
    'R1': (0.18, 0.82), 'R2': (0.82, 0.82),
    'R3': (0.82, 0.18), 'R4': (0.18, 0.18),
}

# ── Layout ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(ns.COL2, 5.2))
gs = gridspec.GridSpec(2, 2, figure=fig,
                       left=0.07, right=0.97,
                       top=0.94, bottom=0.07,
                       wspace=0.38, hspace=0.42,
                       height_ratios=[2.2, 2.0])

axA = fig.add_subplot(gs[0, 0])   # network diagram
axB = fig.add_subplot(gs[0, 1])   # merged PSD
axC = fig.add_subplot(gs[1, :])   # merged time series — full width

ns.label_panel(axA, 'A', x=-0.08, y=1.02)
ns.label_panel(axB, 'B', x=-0.10, y=1.02)
ns.label_panel(axC, 'C', x=-0.04, y=1.02)

# ══════════════════════════════════════════════════════════════════════════════
# Panel A – Network diagram
# ══════════════════════════════════════════════════════════════════════════════
axA.set_xlim(0, 1)
axA.set_ylim(-0.06, 1.06)
axA.set_aspect('equal')
axA.axis('off')

node_r    = 0.10
node_fill = ns.COLORS['pale_grey']
edge_col  = ns.COLORS['dark_grey']

for name, (cx, cy) in node_pos.items():
    circle = mpatches.Circle((cx, cy), node_r,
                             color=node_fill, ec=edge_col,
                             linewidth=0.9, zorder=3,
                             transform=axA.transData, clip_on=False)
    axA.add_patch(circle)
    axA.text(cx, cy, name, ha='center', va='center',
             fontsize=8.4, fontweight='bold', zorder=4)

def _endpoints(src, dst, r=0.10, gap=0.025):
    sx, sy = node_pos[src]
    dx, dy = node_pos[dst]
    ang = np.arctan2(dy - sy, dx - sx)
    return ((sx + (r + gap) * np.cos(ang), sy + (r + gap) * np.sin(ang)),
            (dx - (r + gap) * np.cos(ang), dy - (r + gap) * np.sin(ang)))

for i, src in enumerate(node_labels):
    for j, dst in enumerate(node_labels):
        if i == j:
            continue
        w = Atrue[i, j]
        if w == 0.0:
            continue
        color = ns.COLORS['vermil'] if w > 0 else ns.COLORS['blue']
        lw = 0.5 + 2.8 * abs(w)
        (x0, y0), (x1, y1) = _endpoints(src, dst)
        rad = 0.22 if Atrue[j, i] != 0 else 0.0
        arr = FancyArrowPatch(
            (x0, y0), (x1, y1),
            arrowstyle='-|>', mutation_scale=7,
            connectionstyle=f'arc3,rad={rad}',
            color=color, linewidth=lw, zorder=2,
        )
        axA.add_patch(arr)
        ang_mid = np.arctan2(y1 - y0, x1 - x0)
        mx = (x0 + x1) / 2 + 0.07 * -np.sin(ang_mid) * np.sign(rad if rad else 1)
        my = (y0 + y1) / 2 + 0.07 *  np.cos(ang_mid) * np.sign(rad if rad else 1)
        axA.text(mx, my, f'{w:.1f}', ha='center', va='center',
                 fontsize=6.2, color=color, zorder=5)

for name in node_labels:
    cx, cy = node_pos[name]
    arc_r  = 0.075
    is_top = name in ('R1', 'R2')
    lcy = (cy + node_r + arc_r - 0.01) if is_top else (cy - node_r - arc_r + 0.01)
    ang = np.linspace(0, 2 * np.pi, 80)
    axA.plot(cx + arc_r * np.cos(ang), lcy + arc_r * np.sin(ang),
             color=ns.COLORS['mid_grey'], linewidth=0.8, zorder=1)
    ty = lcy + arc_r + 0.04 if is_top else lcy - arc_r - 0.04
    axA.text(cx, ty, '−0.5', ha='center', va='center',
             fontsize=5.4, color=ns.COLORS['mid_grey'])

leg_A = [
    Line2D([0], [0], color=ns.COLORS['vermil'], lw=1.5, label='Positive ($A_{ij}>0$)'),
    Line2D([0], [0], color=ns.COLORS['blue'],   lw=1.5, label='Negative ($A_{ij}<0$)'),
]
axA.legend(handles=leg_A, loc='lower center', fontsize=7.2,
           frameon=False, handlelength=1.5, bbox_to_anchor=(0.5, -0.04))
axA.set_title('Ground-truth network ($A_{\\rm true}$)', fontsize=8.4, pad=3)

# ══════════════════════════════════════════════════════════════════════════════
# Panel B – Merged PSD (neural solid + Ca²⁺ dashed darker)
# ══════════════════════════════════════════════════════════════════════════════
ns.despine(axB)
for i in range(N):
    axB.plot(_freqs, _psd_x[i], color=colors4[i],  linewidth=0.9, linestyle='-')
    axB.plot(_freqs, _psd_y[i], color=dark4[i],    linewidth=0.9, linestyle='--')

axB.set_xlim(_freqs.min(), _freqs.max())
axB.set_yscale('log')
axB.set_xlabel('Frequency (Hz)', fontsize=7.2)
axB.set_ylabel('Power (a.u.)',   fontsize=7.2)
axB.set_title('Power spectral density', fontsize=8.4, pad=3)
axB.tick_params(labelsize=5.5)

# Legend: nodes by color + line style for signal type
leg_nodes = [Line2D([0], [0], color=colors4[i], lw=1.2, label=node_labels[i])
             for i in range(N)]
leg_style = [
    Line2D([0], [0], color='k', lw=1.0, linestyle='-',  label='Neural state'),
    Line2D([0], [0], color='k', lw=1.0, linestyle='--', label=r'Ca$^{2+}$'),
]
axB.legend(handles=leg_nodes + leg_style, frameon=False, fontsize=6.0,
           ncol=2, loc='upper right', handlelength=1.4, handletextpad=0.4,
           columnspacing=0.8)

# ══════════════════════════════════════════════════════════════════════════════
# Panel C – Merged time series (neural lighter + Ca²⁺ darker, overlaid)
# ══════════════════════════════════════════════════════════════════════════════
ns.despine(axC)
T_show = int(300 / TR)
t_ax   = np.arange(T_show) * TR
offset = 4.5

for i in range(N):
    base = i * offset
    axC.plot(t_ax, x_z[:T_show, i] + base,
             color=colors4[i], linewidth=0.7, alpha=0.55, zorder=2)
    axC.plot(t_ax, y_z[:T_show, i] + base,
             color=dark4[i],   linewidth=0.8, zorder=3)

axC.set_xlabel('Time (s)', fontsize=7.2)
axC.set_yticks([i * offset for i in range(N)])
axC.set_yticklabels(node_labels, fontsize=8.4)
axC.set_xlim(0, T_show * TR)
axC.set_title(r'Latent neural state (lighter) and Ca$^{2\!+}$ observation (darker)', fontsize=8.4, pad=3)
axC.tick_params(axis='x', labelsize=5.5)
axC.tick_params(axis='y', length=0)

leg_C = [
    Line2D([0], [0], color=ns.COLORS['dark_grey'], lw=1.0, alpha=0.55,
           label='Neural state (latent)'),
    Line2D([0], [0], color=darken(ns.COLORS['dark_grey']), lw=1.0,
           label=r'Ca$^{2+}$ observation'),
]
axC.legend(handles=leg_C, frameon=False, fontsize=6.6,
           loc='lower right', bbox_to_anchor=(1.0, 1.01), handlelength=1.4)

# ── Save ──────────────────────────────────────────────────────────────────────
# Shift C right so y-tick labels start at A's left edge; right edge matches B
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
fig_w_px = fig.get_window_extent(renderer).width
bbA = axA.get_window_extent(renderer)
bbB = axB.get_window_extent(renderer)
bbC = axC.get_window_extent(renderer)
yl_bboxes = [t.get_window_extent(renderer) for t in axC.get_yticklabels() if t.get_visible()]
yl_x0_px = min(b.x0 for b in yl_bboxes) if yl_bboxes else bbC.x0
yl_w  = (bbC.x0 - yl_x0_px) / fig_w_px
a_x0  = bbA.x0 / fig_w_px
b_x1  = bbB.x1 / fig_w_px
new_x0 = a_x0 + yl_w
pos_C = axC.get_position()
axC.set_position([new_x0, pos_C.y0, b_x1 - new_x0, pos_C.height])

out_path = os.path.join(OUT_DIR, 'Fig1.png')
ns.save_source_data(out_path, {
    'node_labels': node_labels,
    'TR': TR,
    'A_true': Atrue,
    'neural_state_zscore': x_z,
    'calcium_signal_zscore': y_z,
    'freqs_Hz': _freqs,
    'neural_PSD': _psd_x,
    'calcium_PSD': _psd_y,
})
fig.savefig(out_path, dpi=600)
fig.savefig(out_path.replace('.png', '.svg'))
print(f'Saved: {out_path}')
