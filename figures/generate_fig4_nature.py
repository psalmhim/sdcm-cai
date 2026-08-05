"""
generate_fig4_nature.py

Figure 4: Single-subject DCM-CSD fit and posterior connectivity (Nature style).

Data source: zebra/subject_16_DCM_pc1active_14n.mat (HDF5/v7.3).

Panels:
  a  Auto-spectra grid (2x7) - empirical (blue) vs predicted (vermil dashed).
  b  Posterior mean A-matrix (14x14) heatmap with black cell borders.
  c  Top 20 off-diagonal couplings (Ep +/- SD), coloured by sign.

Output: figures/Fig4.png
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import nature_style as ns

ns.apply()

MAT = os.path.join(os.path.dirname(__file__),
                   'zebra', 'subject_16_DCM_pc1active_14n.mat')
OUT = os.path.join(os.path.dirname(__file__),
                   '..', '..', '..', 'Apps', 'Overleaf',
                   'NIMG-DCM-Ca-CSD', 'figures', 'Fig4.png')

NODE_NAMES = ['lTeO_1', 'lTeO_2', 'lTh', 'lP', 'lPT', 'lHb', 'lpRF',
              'rTeO_1', 'rTeO_2', 'rTh', 'rP', 'rPT', 'rHb', 'rpRF']
N = 14


with h5py.File(MAT, 'r') as f:
    A = np.array(f['DCM_est/Ep/A']).T

    def to_complex(arr):
        arr = np.asarray(arr)
        if arr.dtype.names and 'real' in arr.dtype.names:
            return arr['real'] + 1j * arr['imag']
        return arr.astype(complex)

    csd = to_complex(f['DCM_est/Y/csd'][()])
    Hc = to_complex(f['DCM_est/Hc'][()])
    Hz = np.array(f['DCM_est/Y/Hz']).flatten()

    Cp_sparse = sp.csc_matrix((f['DCM_est/Cp/data'][()],
                               f['DCM_est/Cp/ir'][()],
                               f['DCM_est/Cp/jc'][()]))
    Cp_full = Cp_sparse.toarray()
    Cp_diag = np.diag(Cp_full @ Cp_full.T)
    A_sd = np.sqrt(Cp_diag[:N * N]).reshape(N, N)

node_names = NODE_NAMES
disp_names = ns.subscript_labels(NODE_NAMES)


def reorient(arr):
    arr = np.squeeze(np.asarray(arr))
    if arr.shape == (N, N, len(Hz)):
        return arr
    if arr.shape == (len(Hz), N, N):
        return np.transpose(arr, (1, 2, 0))
    freq_ax = [i for i, s in enumerate(arr.shape) if s == len(Hz)][0]
    return np.moveaxis(arr, freq_ax, -1)


csd = reorient(csd)
Hc = reorient(Hc)

emp_auto = np.real(csd[np.arange(N), np.arange(N), :])
pred_auto = np.real(Hc[np.arange(N), np.arange(N), :])


fig = plt.figure(figsize=(ns.COL2, 5.8))

gs_top = gridspec.GridSpec(2, 7, figure=fig,
                           left=0.07, right=0.98,
                           top=0.96, bottom=0.56,
                           wspace=0.42, hspace=0.68)
gs_bot = gridspec.GridSpec(1, 2, figure=fig,
                           left=0.07, right=0.98,
                           top=0.47, bottom=0.06,
                           wspace=0.38)


emp_min = float(emp_auto.min())
emp_max = float(emp_auto.max())
margin = 0.05 * (emp_max - emp_min if emp_max > emp_min else 1.0)
ylim = (emp_min - margin, emp_max + margin)

spec_axes = []
for idx in range(N):
    r, c = divmod(idx, 7)
    ax = fig.add_subplot(gs_top[r, c])
    ax.plot(Hz, emp_auto[idx], color=ns.COLORS['blue'],
            linewidth=0.75, label='Empirical')
    ax.plot(Hz, pred_auto[idx], color=ns.COLORS['vermil'],
            linewidth=0.75, linestyle='--', label='Predicted')
    ax.set_title(disp_names[idx], fontsize=6.7, fontweight='bold', pad=2)
    ax.set_ylim(ylim)
    ax.tick_params(labelsize=4.4, length=1.5, width=0.4)
    ns.despine(ax)
    if r == 1:
        ax.set_xlabel('Hz', fontsize=6.0, labelpad=1.5)
    if c == 0:
        ax.set_ylabel('Power', fontsize=6.0, labelpad=1.5)
    spec_axes.append(ax)

spec_axes[0].legend(frameon=False, fontsize=5.8, handlelength=1.2,
                    labelspacing=0.15, loc='upper right',
                    bbox_to_anchor=(6.8, 1.4))
ns.label_panel(spec_axes[0], "A", x=-0.55, y=1.30)


axB = fig.add_subplot(gs_bot[0, 0])
cmap_div = ns.diverging_cmap()
eye = np.eye(N, dtype=bool)
A_disp = A.copy()
A_disp[np.diag_indices(N)] = -0.5 * np.exp(np.diag(A))
clim = max(np.abs(A_disp[~eye]).max(), 1e-8)
im = axB.imshow(A_disp, cmap=cmap_div, vmin=-clim, vmax=clim,
                aspect='equal', zorder=1)
axB.set_box_aspect(1)

for k in np.arange(-0.5, N + 0.5, 1.0):
    axB.axhline(k, color='black', linewidth=0.8, zorder=3)
    axB.axvline(k, color='black', linewidth=0.8, zorder=3)
axB.set_xlim(-0.5, N - 0.5)
axB.set_ylim(N - 0.5, -0.5)

axB.set_xticks(np.arange(N))
axB.set_yticks(np.arange(N))
axB.set_xticklabels(disp_names, rotation=90, fontsize=6.0)
axB.set_yticklabels(disp_names, fontsize=6.0)
axB.tick_params(length=0)
axB.set_xlabel('Source', fontsize=7.2)
axB.set_ylabel('Target', fontsize=7.2)
axB.set_title('Posterior mean A-matrix', fontsize=8.4)

cb = plt.colorbar(im, ax=axB, shrink=0.75, pad=0.02)
cb.set_label('Ep(A) (Hz)', fontsize=7.2)
cb.ax.tick_params(labelsize=5.0, length=2)
cb.outline.set_linewidth(1.0)
for sp in axB.spines.values():
    sp.set_visible(True)
    sp.set_linewidth(1.2)
ns.label_panel(axB, "B", x=-0.25, y=1.04)


axC = fig.add_subplot(gs_bot[0, 1])

ti, si = np.where(~eye)
vals = A[ti, si]
sds = A_sd[ti, si]
labels = [f'{disp_names[s]}→{disp_names[t]}' for t, s in zip(ti, si)]

order = np.argsort(np.abs(vals))[::-1]
TOP = 20
sel = order[:TOP]

y = np.arange(len(sel))[::-1]
v_sel = vals[sel]
sd_sel = sds[sel]
lab_sel = [labels[i] for i in sel]
bar_colors = [ns.COLORS['vermil'] if v > 0 else ns.COLORS['blue']
              for v in v_sel]

axC.barh(y, v_sel, color=bar_colors, height=0.68,
         edgecolor='none', zorder=2)
axC.errorbar(v_sel, y, xerr=sd_sel, fmt='none',
             ecolor=ns.COLORS['dark_grey'], elinewidth=0.6,
             capsize=1.5, zorder=3)
axC.axvline(0, color='black', linewidth=0.5, zorder=1)
axC.set_yticks(y)
axC.set_yticklabels(lab_sel, fontsize=6.0)
axC.set_xlabel('Ep(A) ± SD (Hz)', fontsize=7.2)
axC.set_title(f'Top {TOP} off-diagonal couplings', fontsize=8.4)
axC.tick_params(length=1.5, width=0.4, labelsize=5.5)

leg = [Patch(facecolor=ns.COLORS['vermil'], label='Positive effect'),
       Patch(facecolor=ns.COLORS['blue'], label='Negative effect')]
axC.legend(handles=leg, frameon=False, fontsize=6.0,
           handlelength=1.0, labelspacing=0.25, loc='lower right')
ns.despine(axC)
ns.grid(axC, axis='x')
ns.label_panel(axC, "C", x=-0.35, y=1.04)


ns.save_source_data(OUT, {
    'node_names': NODE_NAMES,
    'Hz': Hz,
    'emp_auto': emp_auto,
    'pred_auto': pred_auto,
    'A_posterior_mean': A,
    'A_posterior_sd': A_sd,
    'top_coupling_labels': lab_sel,
    'top_coupling_Ep': v_sel,
    'top_coupling_SD': sd_sel,
})

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=600)
fig.savefig(OUT.replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {OUT}')

from PIL import Image
img = Image.open(OUT)
size_kb = os.path.getsize(OUT) / 1024.0
print(f'Dimensions: {img.size[0]} x {img.size[1]} px')
print(f'File size: {size_kb:.1f} KB')
