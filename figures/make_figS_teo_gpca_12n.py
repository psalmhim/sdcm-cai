"""
make_figS_teo_gpca_12n.py

Supplementary figure: TeO superficial-layer cell selection via global spatial PCA
(12-node model). Reproduces FigS_teo_gPCA_12n.png.

Panels:
  Top rows (L/R TeO x 3 projections):
    Scatter of ALL subjects pooled, cells colored superficial (blue) vs deep (red).
    Projections: ML(Y)-DV(Z),  AP(X)-DV(Z),  AP(X)-ML(Y)
  Top-right: PCA score distribution (pooled L and R across subjects)
  Bottom-left: Per-subject active-cell counts in superficial vs deep
  Bottom-right: Pipeline description text box

Global PCA parameters from teo_global_pca.mat:
  g_dir  = [-0.458, +0.744, +0.487]   (pooled across subjects x hemispheres)
  g_mean = [0.607, 0.587, 0.576]

DV convention: superficial (dorsal) corresponds to smaller DV (Z) values.
High PCA score = deep; low score = superficial.  A per-ROI DV-check flips the
label if needed (should rarely trigger given the fixed g_dir).

Usage (run on server with raw data):
  python3 make_figS_teo_gpca_12n.py
"""

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
RAW_DIR  = Path('/remotenas2/remotedata/share/zebrafish/original_raw_data')
DATA_DIR = Path.home() / 'Dropbox/matlabwork/mnet0.92/dcmcai/zebra'
OUT_DIR  = Path.home() / 'Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures'

SUBJECTS = [12, 13, 14, 15, 16, 17, 18]
TEO_SLOT = {'L': 30, 'R': 66}   # 1-based Kunst-atlas slots

# ── Global PCA parameters ───────────────────────────────────────────────────────
pca = sio.loadmat(str(DATA_DIR / 'teo_global_pca.mat'))
g_dir  = np.array(pca['g_dir']).ravel()     # [3]
g_mean = np.array(pca['g_mean']).ravel()    # [3]
evr    = float(np.array(pca['evr']).ravel()[0])

# ── Colours ────────────────────────────────────────────────────────────────────
C_SUPERF = '#2166AC'   # blue  — superficial, selected
C_DEEP   = '#D6604D'   # red   — deep, excluded

# ── Projections to show ─────────────────────────────────────────────────────────
# (xlabel, ylabel, xi, yi, invert_y)
# invert_y=True when DV is on y-axis so dorsal (small DV) appears at top
PROJS = [
    ('ML (Y)', 'DV (Z)', 1, 2, True),
    ('AP (X)', 'DV (Z)', 0, 2, True),
    ('AP (X)', 'ML (Y)', 0, 1, False),
]


def load_teo_cells(sid, side):
    """
    Load raw xyz and time-series for one subject × hemisphere.
    Returns (xyz_raw [N×3], X [N×T]) or (None, None) if unavailable.
    """
    fpath = RAW_DIR / f'subject_{sid}_data.mat'
    if not fpath.exists():
        print(f'  S{sid}: raw file not found'); return None, None
    try:
        d = sio.loadmat(str(fpath), squeeze_me=True)
    except Exception as e:
        print(f'  S{sid}: load error: {e}'); return None, None

    s_data     = d['s_data'].astype(float)    # [N_cells x T]
    result_idx = d['result_idx']
    cell_xyz   = d['CellXYZ'].astype(float)  # [N_cells x 3]

    slot = TEO_SLOT[side]
    idx  = result_idx[slot - 1].astype(int) - 1   # 0-based

    if len(idx) < 20:
        print(f'  S{sid} {side}TeO: too few cells ({len(idx)})'); return None, None

    xyz = cell_xyz[idx, :]
    X   = s_data[idx, :]
    ok  = np.all(np.isfinite(X), axis=1)
    return xyz[ok, :], X[ok, :]


def gpca_split(xyz_raw, mirror_lr=False):
    """
    Apply global-PCA superficial/deep split.
    Returns (scores, lbl_deep [bool], xyz_n) for the given TeO side.
    """
    lo  = xyz_raw.min(axis=0)
    hi  = xyz_raw.max(axis=0)
    rng = np.maximum(hi - lo, 1e-9)
    xyz_n = (xyz_raw - lo) / rng

    if mirror_lr:
        xyz_n[:, 1] = 1.0 - xyz_n[:, 1]

    scores = (xyz_n - g_mean) @ g_dir

    med = np.median(scores)
    lbl_deep = scores > med   # high score = deep (initial assignment)

    # DV convention check: superficial must have smaller DV (axis 2)
    dv_s = xyz_n[~lbl_deep, 2].mean() if (~lbl_deep).sum() > 0 else 0.0
    dv_d = xyz_n[ lbl_deep, 2].mean() if ( lbl_deep).sum() > 0 else 1.0
    if dv_s > dv_d:
        lbl_deep = ~lbl_deep

    return scores, lbl_deep, xyz_n


# ── Collect data ───────────────────────────────────────────────────────────────
print('Loading subjects...')
pool = {'L': {'xyz': [], 'scores': [], 'lbl_deep': []},
        'R': {'xyz': [], 'scores': [], 'lbl_deep': []}}

# per-subject cell counts for bar chart
counts = {sid: {'Ls': 0, 'Ld': 0, 'Rs': 0, 'Rd': 0} for sid in SUBJECTS}

for sid in SUBJECTS:
    for side in ['L', 'R']:
        xyz_raw, X = load_teo_cells(sid, side)
        if xyz_raw is None:
            continue

        # active cells only (variance > within-TeO median)
        cell_var = X.var(axis=1)
        act_mask = cell_var > np.median(cell_var)
        xyz_act  = xyz_raw[act_mask]

        scores, lbl_deep, xyz_n = gpca_split(xyz_act, mirror_lr=(side == 'R'))

        pool[side]['xyz'].append(xyz_act)
        pool[side]['scores'].append(scores)
        pool[side]['lbl_deep'].append(lbl_deep)

        ns = (~lbl_deep).sum()
        nd =  lbl_deep.sum()
        print(f'  S{sid} {side}TeO: active={len(xyz_act)},  superficial={ns},  deep={nd}')
        if side == 'L':
            counts[sid]['Ls'] = ns; counts[sid]['Ld'] = nd
        else:
            counts[sid]['Rs'] = ns; counts[sid]['Rd'] = nd

# Concatenate pools
for side in ['L', 'R']:
    pool[side]['xyz']      = np.concatenate(pool[side]['xyz'],      axis=0) if pool[side]['xyz']      else np.zeros((0,3))
    pool[side]['scores']   = np.concatenate(pool[side]['scores'],   axis=0) if pool[side]['scores']   else np.zeros(0)
    pool[side]['lbl_deep'] = np.concatenate(pool[side]['lbl_deep'], axis=0) if pool[side]['lbl_deep'] else np.zeros(0, dtype=bool)

# ── Figure layout ──────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7,
    'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
    'axes.linewidth': 0.6, 'pdf.fonttype': 42,
})

fig = plt.figure(figsize=(14, 9))

# Grid: top 2 rows = scatter (L/R × 3projs + hist), bottom row = bar + text
from matplotlib.gridspec import GridSpec
gs = GridSpec(3, 5,
              left=0.06, right=0.98, top=0.93, bottom=0.07,
              hspace=0.45, wspace=0.32,
              height_ratios=[1, 1, 0.85])

fig.suptitle('TeO superficial-layer cell selection — global spatial PCA (12-node model)\n'
             'Blue = superficial (TeO$_s$, selected for 12-node signal)  '
             'Red = deep (TeO$_{deep}$, excluded)',
             fontsize=10, y=0.98)

# scatter axes (2 rows × 3 cols)
scat_ax = [[fig.add_subplot(gs[row, col]) for col in range(3)] for row in range(2)]
# histogram axis (right, spans both scatter rows)
hist_ax = fig.add_subplot(gs[0:2, 3])
# bar chart axis (bottom-left, spans 2 cols)
bar_ax  = fig.add_subplot(gs[2, 0:3])
# text box axis (bottom-right)
txt_ax  = fig.add_subplot(gs[2, 3:5])

SIDE_LABELS = ['Left TeO', 'Right TeO']
ALPHA = 0.15
S_SIZE = 0.8

for si, side in enumerate(['L', 'R']):
    xyz    = pool[side]['xyz']
    lbl_d  = pool[side]['lbl_deep']
    mS     = ~lbl_d
    mD     =  lbl_d

    for pi, (xlabel, ylabel, xi, yi, inv_y) in enumerate(PROJS):
        ax = scat_ax[si][pi]
        ax.scatter(xyz[mD, xi], xyz[mD, yi], s=S_SIZE, c=C_DEEP,
                   alpha=ALPHA, lw=0, rasterized=True, zorder=1)
        ax.scatter(xyz[mS, xi], xyz[mS, yi], s=S_SIZE, c=C_SUPERF,
                   alpha=ALPHA, lw=0, rasterized=True, zorder=2)
        if inv_y:
            ax.invert_yaxis()
        ax.set_xlabel(xlabel, fontsize=6.5, labelpad=1)
        ax.set_ylabel(ylabel, fontsize=6.5, labelpad=1)
        if pi == 0:
            ax.set_ylabel(f'{SIDE_LABELS[si]}\n{ylabel}', fontsize=6.5, labelpad=1)
        ax.tick_params(length=2)

# ── Histogram of PCA scores ───────────────────────────────────────────────────
for side, ls, lbl in [('L', '-', 'L superficial'), ('L', '--', 'L deep'),
                       ('R', '-',  'R superficial'), ('R', '--', 'R deep')]:
    scores = pool[side]['scores']
    lbl_d  = pool[side]['lbl_deep']
    if side == 'L':
        mask = ~lbl_d if 'superficial' in lbl else lbl_d
        col  = C_SUPERF if 'superficial' in lbl else C_DEEP
        hist_ax.hist(scores[mask], bins=60, density=True,
                     alpha=0.4, color=col, histtype='stepfilled',
                     label=lbl, linewidth=0)
    else:
        mask = ~lbl_d if 'superficial' in lbl else lbl_d
        col  = C_SUPERF if 'superficial' in lbl else C_DEEP
        hist_ax.hist(scores[mask], bins=60, density=True,
                     alpha=0.3, color=col, histtype='step',
                     linestyle=ls, linewidth=1.0, label=lbl)

hist_ax.axvline(0, color='k', ls=':', lw=0.8)
hist_ax.set_xlabel('Global PCA score\n(xyz − x̄) · g$_{dir}$', fontsize=7)
hist_ax.set_ylabel('Density', fontsize=7)
hist_ax.set_title(f'PCA score distribution\n(superficial vs. deep)',
                  fontsize=7)
hist_ax.text(0.05, 0.97,
             f'g$_{{dir}}$ = [{g_dir[0]:.3f}, {g_dir[1]:.3f}, {g_dir[2]:.3f}]\n'
             f'EVR = {evr*100:.1f}%',
             transform=hist_ax.transAxes, va='top', ha='left',
             fontsize=6, bbox=dict(fc='white', alpha=0.7, pad=2, ec='none'))
hist_ax.legend(fontsize=5.5, loc='upper right', framealpha=0.7)

# ── Bar chart: per-subject cell counts ─────────────────────────────────────────
x_pos  = np.arange(len(SUBJECTS))
width  = 0.2
sids   = list(counts.keys())

bar_ax.bar(x_pos - 1.5*width, [counts[s]['Ls'] for s in sids], width,
           color=C_SUPERF, label='L superficial (A)')
bar_ax.bar(x_pos - 0.5*width, [counts[s]['Ld'] for s in sids], width,
           color='#92C5DE', label='L deep (B)')
bar_ax.bar(x_pos + 0.5*width, [counts[s]['Rs'] for s in sids], width,
           color=C_DEEP, label='R superficial (A)')
bar_ax.bar(x_pos + 1.5*width, [counts[s]['Rd'] for s in sids], width,
           color='#F4A582', label='R deep (B)')
bar_ax.axhline(0, color='k', lw=0.5)
bar_ax.set_xticks(x_pos)
bar_ax.set_xticklabels([f'S{s}' for s in sids], fontsize=7)
bar_ax.set_ylabel('Number of cells', fontsize=7)
bar_ax.set_title('Per-subject active cell counts in consensus cores', fontsize=7)
bar_ax.legend(fontsize=5.5, ncol=4, loc='upper right', framealpha=0.7)

# ── Text box: pipeline description ─────────────────────────────────────────────
txt_ax.axis('off')
pipeline_text = (
    'TeO signal extraction pipeline (12-node)\n\n'
    u'①  Load atlas-registered cell coordinates\n'
    '    (CellXYZ) + calcium time series (s_data)\n\n'
    u'②  Select active cells\n'
    '    (temporal variance > within-ROI median)\n\n'
    u'③  Project onto global PCA axis g₍ᵉᴼᵣ₎\n'
    '    (fit once on all subjects × hemispheres)\n'
    '    Score s = (xyz − x̅) · gₚᵈᵃ\n'
    '    Split at per-ROI-score median\n'
    '    TeOₛ = low-score half (superficial, smaller DV)\n'
    '    TeO₞ = high score half (excluded)\n\n'
    u'④  Per-cell standardisation\n'
    '    detrend → mean-centre → unit-variance\n\n'
    u'⑤  Population signal = mean of standardised cells\n'
    '    (sign-aligned to ROI mean)\n'
)
txt_ax.text(0.02, 0.98, pipeline_text,
            transform=txt_ax.transAxes, va='top', ha='left',
            fontsize=6.5, family='monospace',
            bbox=dict(fc='#F8F8F8', ec='#CCCCCC', pad=6, lw=0.6))

# ── Panel labels ───────────────────────────────────────────────────────────────
for ax, lbl in zip([scat_ax[0][0], scat_ax[1][0], hist_ax, bar_ax],
                   ['A', 'B', 'C', 'D']):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.01, pos.y1 + 0.005, lbl,
             fontsize=10, fontweight='bold', va='bottom')

# ── Legend patches at figure level ─────────────────────────────────────────────
fig.legend(handles=[mpatches.Patch(color=C_SUPERF, label='Superficial (TeO$_s$, selected)'),
                    mpatches.Patch(color=C_DEEP,   label='Deep (TeO$_{deep}$, excluded)')],
           loc='lower center', ncol=2, fontsize=8, bbox_to_anchor=(0.3, 0.002),
           framealpha=0.9)

# ── Save ───────────────────────────────────────────────────────────────────────
OUT_DIR.mkdir(parents=True, exist_ok=True)
out = OUT_DIR / 'FigS_teo_gPCA_12n.png'
fig.savefig(str(out), dpi=200, bbox_inches='tight')
print(f'\nSaved: {out}')
plt.close()
