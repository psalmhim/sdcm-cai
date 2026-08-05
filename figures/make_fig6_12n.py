"""
make_fig6_12n.py
Group-level effective connectivity, 12-node RegLin-Mean primary specification,
uniform regularized prior (p_C=1/16), post self-connection-prior fix (2026-08-04).

Layout: 1 row x 2 columns
  A: PEB Ep[A], thresholded by Gaussian Pp > 0.975
  B: BMR/BMA Ep[A], confirmed at posterior inclusion probability Pp > 0.975 (=1.0)

Conventions (matching manuscript Fig. 6 caption):
  Off-diagonal A(i,j): directed influence from source j to target i.
  Rows = target, columns = source.
  Diagonal displayed in physical self-inhibition space: -0.5*exp(A_ii).
Output: figures/Fig7_p16.png / .pdf / .svg
"""

import os
import scipy.io
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.expanduser('~/Dropbox/matlabwork/mnet0.92/dcmcai/zebra/')
FIGS = os.path.expanduser('~/Dropbox/Apps/Overleaf/NIMG-DCM-Ca-CSD/figures/')
os.makedirs(FIGS, exist_ok=True)

mpl.rcParams.update({
    'font.family':    'Helvetica',
    'font.size':       7.2,
    'axes.labelsize':  8.4,
    'axes.titlesize':  7.8,
    'xtick.labelsize': 6.0,
    'ytick.labelsize': 6.0,
    'axes.linewidth':  0.6,
    'xtick.major.width': 0.6, 'xtick.major.size': 2.5,
    'ytick.major.width': 0.6, 'ytick.major.size': 2.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

n = 12
short = ['lTeO', 'lTh', 'lP', 'lPT', 'lHb', 'lpRF',
         'rTeO', 'rTh', 'rP', 'rPT', 'rHb', 'rpRF']
offdiag = ~np.eye(n, dtype=bool)
diag_idx = np.arange(n)
THR = 0.975


def ep_disp(Ep):
    E = Ep.copy()
    E[diag_idx, diag_idx] = -0.5 * np.exp(Ep[diag_idx, diag_idx])
    return E


def cell_borders(ax, lw=0.7):
    for k in np.arange(-0.5, n + 0.5, 1.0):
        ax.axhline(k, color='black', lw=lw, zorder=3)
        ax.axvline(k, color='black', lw=lw, zorder=3)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(1.0)


def mark_sig(ax, SIG, Ep_d):
    for r in range(n):
        for c in range(n):
            if r == c:
                continue
            if SIG[r, c]:
                color = 'white' if Ep_d[r, c] < -0.2 else 'black'
                ax.add_patch(plt.Rectangle((c - .5, r - .5), 1, 1,
                    lw=0.6, edgecolor=color, facecolor='none', zorder=3))
                ax.text(c, r + 0.2, '*', ha='center', va='center',
                        fontsize=5.4, color=color, fontweight='bold')


def draw_panel(fig, ax, Ep_d, SIG, title, n_sig, ylabel=True, xlabel=True):
    vmax = np.ceil(abs(Ep_d[offdiag]).max() * 10) / 10
    if vmax < 0.1:
        vmax = 0.1
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(Ep_d, cmap='RdBu_r', norm=norm, aspect='equal',
                    interpolation='nearest')
    mark_sig(ax, SIG, Ep_d)
    ax.axhline(5.5, color='black', lw=0.9)
    ax.axvline(5.5, color='black', lw=0.9)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    if xlabel:
        ax.set_xticklabels(short, rotation=90, fontsize=5.8)
        ax.set_xlabel('Source', labelpad=1, fontsize=7.6)
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', length=0, bottom=False, labelbottom=False)
    if ylabel:
        ax.set_yticklabels(short, fontsize=5.8)
        ax.set_ylabel('Target', labelpad=1, fontsize=7.6)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis='y', length=0, left=False, labelleft=False)
    ax.set_title(f'{title}\n($n_{{sig}}$ = {n_sig}, $P_p$ > {THR})',
                 fontsize=6.8, fontweight='bold', pad=3)
    cell_borders(ax)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, shrink=0.85)
    cb.set_label('$E_p$ (a.u.)', fontsize=5.4, labelpad=1)
    cb.ax.tick_params(labelsize=4.0, width=0.5, length=1.5)


d = scipy.io.loadmat(DATA + 'fig6_arrays_12n.mat', squeeze_me=True)
Ep_peb = d['Ep_peb'].real
Pp_peb = d['Pp_peb'].real
Ep_bma = d['Ep_bma'].real
Pp_bma = d['Pp_bma'].real

SIG_peb = offdiag & (Pp_peb > THR)
SIG_bma = offdiag & (Pp_bma > THR)
n_peb = int(SIG_peb.sum())
n_bma = int(SIG_bma.sum())
print(f'PEB sig: {n_peb}   BMA sig: {n_bma}')

fig, axes = plt.subplots(1, 2, figsize=(5.5, 3.4))
fig.subplots_adjust(left=0.12, right=0.97, bottom=0.14, top=0.86, wspace=0.28)

for idx, ax in enumerate(axes):
    lbl = 'AB'[idx]
    ax.text(-0.06, 1.02, lbl, transform=ax.transAxes,
            fontsize=10.8, fontweight='bold', va='bottom', ha='left')

draw_panel(fig, axes[0], ep_disp(Ep_peb), SIG_peb,
           'PEB $E_p[A]$ (RegLin-Mean, $p_C=1/16$)',
           n_sig=n_peb, ylabel=True, xlabel=True)
draw_panel(fig, axes[1], ep_disp(Ep_bma), SIG_bma,
           'BMR/BMA $E_p[A]$',
           n_sig=n_bma, ylabel=False, xlabel=True)

fig.suptitle('Group-level effective connectivity (RegLin-Mean, uniform prior $p_C=1/16$)',
             fontsize=8.4, fontweight='bold', y=0.97)

for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(FIGS, f'Fig7_p16.{ext}')
    fig.savefig(path, dpi=600, bbox_inches='tight')
    print(f'Saved {path}')
plt.close()
