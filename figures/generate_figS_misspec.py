"""
generate_figS_misspec.py

Regenerate FigS_misspec.png (Supplementary Figure S1) showing only the
nonzero-connection scatter (red). Blue "All 132 off-diag" scatter removed.

Output: figures/FigS_misspec.png (Overleaf)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.io
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nature_style as ns

ns.apply()

HERE    = os.path.dirname(os.path.abspath(__file__))
MAT     = os.path.join(HERE, "zebra", "misspecified_model_results.mat")
OUT_DIR = os.path.join(HERE, "..", "..", "..", "Apps", "Overleaf", "NIMG-DCM-Ca-CSD", "figures")
OUT     = os.path.join(OUT_DIR, "FigS_misspec.png")

mat      = scipy.io.loadmat(MAT, squeeze_me=True, struct_as_record=False)
A_true   = mat["A_true"]          # (12,12) — 12-node synthetic network
A_est    = mat["DCM_est"].Ep.A    # (12,12)
r_nz     = float(mat["r_nonzero"])
sign_nz  = float(mat["sign_nonzero"])
n_nz     = int(mat["n_nonzero"])
n        = A_true.shape[0]

print(f'=== Misspecification simulation (12-node, linear SDE → sDCM-CaI) ===')
print(f'  r (nonzero) = {r_nz:.3f}  (Paper: 0.932)')
print(f'  sign acc    = {sign_nz*100:.1f}%  (Paper: 96.2%)')

# Nonzero off-diagonal mask (from A_true)
offdiag = ~np.eye(n, dtype=bool)
nonzero = offdiag & (A_true != 0)

v_true_nz = A_true[nonzero]
v_rec_nz  = A_est[nonzero]

# Colour by true sign
col_nz = np.where(v_true_nz > 0, ns.COLORS["vermil"], ns.COLORS["blue"])

fig, ax = plt.subplots(figsize=(3.5, 3.5))

lim = (-0.45, 0.45)
ax.plot(lim, lim, "k--", linewidth=0.8, zorder=1, alpha=0.6, label="Identity")

ax.scatter(v_true_nz, v_rec_nz,
           s=40, c=col_nz,
           edgecolors=ns.COLORS["ink"], linewidths=0.4, zorder=4)

ax.set_xlim(*lim)
ax.set_ylim(*lim)
ax.set_aspect("equal")
ax.set_xlabel(r"True $A_{ij}$ (linear SDE)", fontsize=7.8)
ax.set_ylabel(r"Recovered $A_{ij}$ (sDCM-CaI, misspecified)", fontsize=7.8)
ax.set_title("Misspecified-model recovery", fontsize=9.0)
ax.tick_params(labelsize=5.5)

ax.text(0.05, 0.94,
        f"Nonzero {n_nz}: $r={r_nz:.3f}$, sign$={sign_nz*100:.1f}\\%$",
        transform=ax.transAxes, fontsize=7.2,
        va="top", color=ns.COLORS["vermil"])

legend_elems = [
    Line2D([0], [0], color="k", ls="--", lw=0.8, label="Identity"),
    Line2D([0], [0], marker="o", color="w",
           mfc=ns.COLORS["vermil"], mec=ns.COLORS["ink"], ms=5, label="Positive effect (true > 0)"),
    Line2D([0], [0], marker="o", color="w",
           mfc=ns.COLORS["blue"],   mec=ns.COLORS["ink"], ms=5, label="Negative effect (true < 0)"),
]
ax.legend(handles=legend_elems, frameon=False, fontsize=6.2,
          loc="lower right", handlelength=1.2, handletextpad=0.4)

ns.despine(ax)
ns.grid(ax)

ns.save_source_data(OUT, {
    'A_true': A_true,
    'A_recovered': A_est,
    'true_nonzero_offdiag': v_true_nz,
    'recovered_nonzero_offdiag': v_rec_nz,
    'r_nonzero': r_nz,
    'sign_nonzero': sign_nz,
    'n_nonzero': n_nz,
})

fig.tight_layout()
os.makedirs(OUT_DIR, exist_ok=True)
fig.savefig(OUT, dpi=600, bbox_inches="tight")
fig.savefig(OUT.replace(".png", ".svg"), bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")
