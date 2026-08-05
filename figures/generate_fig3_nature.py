"""
Generate Figure 3 (saved as figures/Fig3.png) - Nature style.

Panel A (left column, stacked vertically):
  A-top:    14×14 structural connectivity (SC) matrix from Kunst atlas
  A-bottom: 14×14 Pearson functional connectivity (FC) matrix
Panel B (right): stacked z-scored PC1 calcium time series.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import h5py
import scipy.io as sio

sys.path.insert(0, os.path.dirname(__file__))
import nature_style as ns

ns.apply()

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "zebra", "subject_16_pc1signals_14n.mat")
SC_MAT = os.path.join(HERE, "zebrafish_sc_mat.mat")
OUT = os.path.join(
    HERE, "..", "..", "..",
    "Apps", "Overleaf", "NIMG-DCM-Ca-CSD", "figures", "Fig3.png",
)

TR = 0.5
NSHOW = None   # None = show all timepoints
OFFSET = 4.5
HEMI_SPLIT = 6.5   # between lpRF (idx 6) and rTeO_s (idx 7)

NODE_NAMES = [
    "lTeO_1", "lTeO_2", "lTh", "lP", "lPT", "lHb", "lpRF",
    "rTeO_1", "rTeO_2", "rTh", "rP", "rPT", "rHb", "rpRF",
]

# 14-node atlas indices (1-based MATLAB → 0-based Python)
SC_IDX = [i - 1 for i in [30, 30, 31, 23, 25, 17, 15, 66, 66, 67, 59, 61, 53, 51]]

REGION_COLORS = {
    'TeO_1': '#D8A33D',
    'TeO_2': '#0000FF',
    'Th':    '#1B9A89',
    'P':     '#C65A2E',
    'PT':    '#CC79A7',
    'Hb':    '#2F6F9F',
    'pRF':   '#7F7F7F',
}


def get_node_color(name):
    region = name[1:]
    return REGION_COLORS.get(region, '#888888')


def load_signals():
    with h5py.File(DATA, "r") as f:
        sig = np.array(f["pc1_active_signals"])
        if sig.shape[0] == len(NODE_NAMES):
            sig = sig.T
    return sig.astype(float)


def load_sc14():
    sc = sio.loadmat(SC_MAT, squeeze_me=True)
    SC72 = sc['sc_mat'].astype(float)
    SC14 = SC72[np.ix_(SC_IDX, SC_IDX)]
    np.fill_diagonal(SC14, 0)   # zero self-connections for display
    return SC14


def _add_matrix_grid(ax, n):
    for k in np.arange(-0.5, n + 0.5, 1.0):
        ax.axhline(k, color='black', linewidth=0.3, zorder=3)
        ax.axvline(k, color='black', linewidth=0.3, zorder=3)
    ax.axhline(HEMI_SPLIT, color='black', linewidth=0.9, zorder=4)
    ax.axvline(HEMI_SPLIT, color='black', linewidth=0.9, zorder=4)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(ns.subscript_labels(NODE_NAMES), rotation=90, fontsize=4.5)
    ax.set_yticklabels(ns.subscript_labels(NODE_NAMES), fontsize=4.5)
    ax.tick_params(length=0)
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel('Source', fontsize=6.5)
    ax.set_ylabel('Target', fontsize=6.5)


def panel_a_sc(ax, SC14):
    n = 14
    sc_log = np.log1p(SC14)
    im = ax.imshow(sc_log, cmap='YlOrRd', aspect='equal', interpolation='nearest')
    _add_matrix_grid(ax, n)
    ax.set_title('Structural connectivity (SC)', fontsize=8.0)
    cb = ax.get_figure().colorbar(im, ax=ax, orientation='horizontal',
                                  fraction=0.046, pad=0.20, shrink=0.85)
    cb.set_label('Fibre count (log)', fontsize=5.0, labelpad=1)
    cb.ax.tick_params(labelsize=4.0, width=0.5, length=2)
    ns.label_panel(ax, 'A', x=-0.22, y=1.06)


def panel_a_fc(ax, sig):
    n = 14
    FC = np.corrcoef(sig[:, :n].T)
    im = ax.imshow(FC, cmap=ns.diverging_cmap(), vmin=-1, vmax=1,
                   aspect='equal', interpolation='nearest', zorder=1)
    _add_matrix_grid(ax, n)
    ax.set_title('Functional connectivity (FC)', fontsize=8.0)
    cb = ax.get_figure().colorbar(im, ax=ax, orientation='horizontal',
                                  fraction=0.046, pad=0.20, shrink=0.85)
    cb.set_label('Pearson r', fontsize=5.0, labelpad=1)
    cb.set_ticks([-1, -0.5, 0, 0.5, 1])
    cb.ax.tick_params(labelsize=4.0, width=0.5, length=2)


def panel_b(ax, sig, names):
    nshow = sig.shape[0] if NSHOW is None else min(NSHOW, sig.shape[0])
    t = np.arange(nshow) * TR
    n = len(names)
    yticks = []

    for i in range(n):
        x = sig[:nshow, i].astype(float)
        sd = np.std(x)
        if sd > 0:
            x = (x - np.mean(x)) / sd
        base = (n - 1 - i) * OFFSET
        color = get_node_color(names[i])
        ax.plot(t, x + base, color=color, linewidth=0.55, zorder=2)
        yticks.append(base)

    ax.set_yticks(yticks)
    ax.set_yticklabels(ns.subscript_labels(names), fontsize=5.4)
    max_t = int(round(nshow * TR))
    ax.set_xlim(0, max_t)
    ax.set_ylim(-OFFSET, n * OFFSET)
    xticks = list(range(0, max_t, 50)) + [max_t]
    ax.set_xticks(xticks)
    ax.set_xlabel("Time (s)", fontsize=7.2)
    ax.set_title("PC1 calcium signals (subj. 16)", fontsize=8.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=5.0)
    ns.label_panel(ax, "B", x=-0.14, y=1.02)


def main():
    sig = load_signals()
    SC14 = load_sc14()
    names = list(NODE_NAMES)
    print(f"Signals: {sig.shape}  SC14: {SC14.shape}")

    fig = plt.figure(figsize=(ns.COL2, 5.5))

    # Panel A: SC + FC matrices stacked in left column
    gs_a = gridspec.GridSpec(2, 1, figure=fig,
                             left=0.07, right=0.38,
                             top=0.96, bottom=0.05,
                             hspace=0.60)
    # Panel B: time series in right column
    gs_b = gridspec.GridSpec(1, 1, figure=fig,
                             left=0.45, right=0.97,
                             top=0.96, bottom=0.05)

    axA_sc = fig.add_subplot(gs_a[0])
    axA_fc = fig.add_subplot(gs_a[1])
    axB    = fig.add_subplot(gs_b[0])

    panel_a_sc(axA_sc, SC14)
    panel_a_fc(axA_fc, sig)
    panel_b(axB, sig, names)

    n = len(names)
    ns.save_source_data(OUT, {
        'node_names': names,
        'TR': TR,
        'pc1_signals': sig[:, :n],
        'SC14': SC14,
        'FC14': np.corrcoef(sig[:, :n].T),
    })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=600, bbox_inches='tight')
    fig.savefig(OUT.replace('.png', '.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {OUT}")

    try:
        from PIL import Image
        with Image.open(OUT) as img:
            w, h = img.size
        print(f"Image: {w}×{h} px  ({os.path.getsize(OUT)/1024:.1f} KB)")
    except Exception as e:
        print(f"PIL check failed: {e}")


if __name__ == "__main__":
    main()
