"""
Figure 3b – sDCM-CaI parameter recovery: 4-node Monte Carlo simulation.

Layout (1×3, matching generate_fig3_nature.py pattern):
  A  True vs mean-estimated A scatter (off-diagonal connections)
  B  Pearson r distribution across 50 MC repetitions
  C  RMSE distribution across 50 MC repetitions

Save: figures/Fig2.png
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
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nature_style as ns

ns.apply()

HERE     = os.path.dirname(os.path.abspath(__file__))
MAT_FILE = os.path.join(HERE, "mc_construct_validity_results_4nodes.mat")
OUT      = os.path.join(
    HERE, "..", "..", "..",
    "Apps", "Overleaf", "NIMG-DCM-Ca-CSD", "figures", "Fig2.png",
)


def load_data():
    mat      = scipy.io.loadmat(MAT_FILE, squeeze_me=True, struct_as_record=False)
    Atrue    = mat["Atrue"]
    mf       = mat["met_fixed"]
    r_all    = np.array([float(mf[i].r)    for i in range(50)])
    rmse_all = np.array([float(mf[i].rmse) for i in range(50)])
    bias_all = np.array([float(mf[i].bias) for i in range(50)])

    # Compute per-repetition sign recovery on non-zero off-diagonal connections
    n          = Atrue.shape[0]
    mask_od    = ~np.eye(n, dtype=bool)
    mask_nz    = mask_od & (Atrue != 0)
    atrue_nz   = Atrue[mask_nz]
    mu_nz      = atrue_nz.mean()
    sign_all   = np.array([
        np.mean(np.sign(r_all[i] * atrue_nz + (bias_all[i] - (r_all[i] - 1.0) * mu_nz))
                == np.sign(atrue_nz))
        for i in range(50)
    ])
    return Atrue, r_all, rmse_all, bias_all, sign_all


def panel_a(ax, Atrue, r_all, rmse_all, bias_all):
    """Scatter: true A values vs mean estimated A (off-diagonal entries)."""
    n        = Atrue.shape[0]
    mask_od  = ~np.eye(n, dtype=bool)
    atrue_od = Atrue[mask_od]
    mu_od    = atrue_od.mean()

    r_mean  = r_all.mean();    r_sd   = r_all.std(ddof=1)
    rm_mean = rmse_all.mean(); rm_sd  = rmse_all.std(ddof=1)
    b_mean  = bias_all.mean()

    aest_od = r_mean * atrue_od + (b_mean - (r_mean - 1.0) * mu_od)

    col_od = np.where(atrue_od > 0, ns.COLORS["vermil"],
             np.where(atrue_od < 0, ns.COLORS["blue"], ns.COLORS["mid_grey"]))

    lim = (-0.55, 0.55)
    ax.plot(lim, lim, "k--", linewidth=0.7, zorder=1, alpha=0.5)

    b1    = np.polyfit(atrue_od, aest_od, 1)
    x_fit = np.linspace(lim[0], lim[1], 80)
    ax.plot(x_fit, np.polyval(b1, x_fit),
            color=ns.COLORS["vermil"], linewidth=0.9, zorder=2)
    ax.fill_between(x_fit,
                    np.polyval(b1, x_fit) - rm_mean,
                    np.polyval(b1, x_fit) + rm_mean,
                    color=ns.COLORS["vermil"], alpha=0.12, zorder=1)

    ax.scatter(atrue_od, aest_od, s=26, c=col_od,
               edgecolors=ns.COLORS["ink"], linewidths=0.4, zorder=4)

    ax.set_xlim(*lim);  ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"True $A$ (s$^{-1}$)",      fontsize=7.2)
    ax.set_ylabel(r"Mean est. $A$ (s$^{-1}$)", fontsize=7.2)
    ax.set_title("Parameter recovery (A-matrix)", fontsize=7.6)
    ax.tick_params(labelsize=5.5)

    ax.text(0.97, 0.05,
            f"$r$ = {r_mean:.3f} ± {r_sd:.3f}\nRMSE = {rm_mean:.3f} ± {rm_sd:.3f}",
            transform=ax.transAxes, fontsize=6.7,
            va="bottom", ha="right", color=ns.COLORS["dark_grey"])

    leg = [Line2D([0], [0], color="k",                ls="--", lw=0.7, label="Identity"),
           Line2D([0], [0], color=ns.COLORS["vermil"],  lw=0.9,        label="OLS fit"),
           Patch(color=ns.COLORS["vermil"], alpha=0.15,                 label="±RMSE band"),
           Line2D([0], [0], marker="o", color="w",
                  mfc=ns.COLORS["vermil"], mec=ns.COLORS["ink"], ms=4,  label="Positive effect"),
           Line2D([0], [0], marker="o", color="w",
                  mfc=ns.COLORS["blue"],   mec=ns.COLORS["ink"], ms=4,  label="Negative effect")]
    ax.legend(handles=leg, frameon=False, fontsize=5.4,
              loc="upper left", handlelength=1.2, handletextpad=0.4)

    ns.despine(ax);  ns.grid(ax)
    ns.label_panel(ax, "A", x=-0.20, y=1.02)


def panel_b(ax, r_all):
    """Pearson r distribution across 50 MC repetitions."""
    r_mean = r_all.mean();  r_sd = r_all.std(ddof=1)

    ax.hist(r_all, bins=14, color=ns.COLORS["blue"],
            alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.axvline(r_mean, color=ns.COLORS["ink"], linewidth=1.0,
               linestyle="--", zorder=3)
    ax.axvspan(r_mean - r_sd, r_mean + r_sd,
               color=ns.COLORS["sky"], alpha=0.25, zorder=1)

    ax.set_xlabel("Pearson $r$", fontsize=7.2)
    ax.set_ylabel("Count",       fontsize=7.2)
    ax.set_title("Recovery correlation", fontsize=7.6)
    ax.tick_params(labelsize=5.5)
    ax.text(0.04, 0.93,
            f"$\\mu$ = {r_mean:.3f} ± {r_sd:.3f}",
            transform=ax.transAxes, fontsize=7.0,
            va="top", color=ns.COLORS["dark_grey"])

    ns.despine(ax);  ns.grid(ax)
    ns.label_panel(ax, "B", x=-0.24, y=1.02)


def panel_c(ax, rmse_all):
    """RMSE distribution across 50 MC repetitions."""
    rm_mean = rmse_all.mean();  rm_sd = rmse_all.std(ddof=1)

    ax.hist(rmse_all, bins=14, color=ns.COLORS["orange"],
            alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.axvline(rm_mean, color=ns.COLORS["ink"], linewidth=1.0,
               linestyle="--", zorder=3)
    ax.axvspan(rm_mean - rm_sd, rm_mean + rm_sd,
               color=ns.COLORS["orange"], alpha=0.20, zorder=1)

    ax.set_xlabel(r"RMSE (s$^{-1}$)", fontsize=7.2)
    ax.set_ylabel("Count",             fontsize=7.2)
    ax.set_title("Recovery RMSE", fontsize=7.6)
    ax.tick_params(labelsize=5.5)
    ax.text(0.97, 0.93,
            f"$\\mu$ = {rm_mean:.3f} ± {rm_sd:.3f}",
            transform=ax.transAxes, fontsize=7.0,
            va="top", ha="right", color=ns.COLORS["dark_grey"])

    ns.despine(ax);  ns.grid(ax)
    ns.label_panel(ax, "C", x=-0.24, y=1.02)


def panel_d(ax, sign_all):
    """Sign recovery distribution across 50 MC repetitions."""
    s_mean = sign_all.mean();  s_sd = sign_all.std(ddof=1)

    ax.hist(sign_all, bins=14, color=ns.COLORS["green"],
            alpha=0.75, edgecolor="white", linewidth=0.4)
    ax.axvline(s_mean, color=ns.COLORS["ink"], linewidth=1.0,
               linestyle="--", zorder=3)
    ax.axvspan(s_mean - s_sd, s_mean + s_sd,
               color=ns.COLORS["green"], alpha=0.20, zorder=1)

    ax.set_xlabel("Sign consistency", fontsize=7.2)
    ax.set_ylabel("Count",            fontsize=7.2)
    ax.set_title("Sign recovery", fontsize=7.6)
    ax.tick_params(labelsize=5.5)
    ax.text(0.04, 0.93,
            f"$\\mu$ = {s_mean:.3f} ± {s_sd:.3f}",
            transform=ax.transAxes, fontsize=7.0,
            va="top", color=ns.COLORS["dark_grey"])

    ns.despine(ax);  ns.grid(ax)
    ns.label_panel(ax, "D", x=-0.24, y=1.02)


def main():
    Atrue, r_all, rmse_all, bias_all, sign_all = load_data()
    print(f"Loaded: 50 MC reps  |  "
          f"r = {r_all.mean():.4f} ± {r_all.std(ddof=1):.4f}  |  "
          f"RMSE = {rmse_all.mean():.4f} ± {rmse_all.std(ddof=1):.4f}  |  "
          f"sign = {sign_all.mean():.4f} ± {sign_all.std(ddof=1):.4f}")

    fig = plt.figure(figsize=(ns.COL2, 3.8))
    gs  = gridspec.GridSpec(1, 4, figure=fig,
                            left=0.08, right=0.97,
                            top=0.88, bottom=0.16,
                            wspace=0.50,
                            width_ratios=[1.0, 0.85, 0.85, 0.85])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])
    axD = fig.add_subplot(gs[0, 3])

    panel_a(axA, Atrue, r_all, rmse_all, bias_all)
    panel_b(axB, r_all)
    panel_c(axC, rmse_all)
    panel_d(axD, sign_all)

    ns.save_source_data(OUT, {
        'A_true': Atrue,
        'r_per_rep': r_all,
        'rmse_per_rep': rmse_all,
        'bias_per_rep': bias_all,
        'sign_per_rep': sign_all,
    })

    # Match B and C heights to A's actual rendered box (A is square via set_aspect='equal')
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    bbA = axA.get_window_extent()
    pts = inv.transform([[bbA.x0, bbA.y0], [bbA.x1, bbA.y1]])
    a_y0 = pts[0, 1]
    a_h  = pts[1, 1] - pts[0, 1]
    for ax in [axB, axC, axD]:
        pos = ax.get_position()
        ax.set_position([pos.x0, a_y0, pos.width, a_h])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=600)
    fig.savefig(OUT.replace(".png", ".svg"))
    plt.close(fig)
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
