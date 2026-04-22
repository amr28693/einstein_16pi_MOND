"""
Generate the headline empirical figure for the MOND-BE paper by Rodriguez (2026).

Three-panel figure:
  (a) D_obs vs theta, with Bose-Einstein prediction curve overlaid
      and the ln(2) crossover marked at (ln 2, 2).
  (b) Residuals in log-space (log10 D_obs - log10 D_pred) vs gbar_log,
      with mean and +/- sigma bands.
  (c) n_BE as a function of theta from data (binned medians), compared
      to the prediction curve, spanning the full range.

Uses Rodriguez_2026_MOND_BE empirical results exactly as reported:
  N = 2693, residual_std = 0.133 dex, crossover at theta = ln 2.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "text.usetex": False,
    "mathtext.fontset": "cm",
})

df = pd.read_csv("bose_einstein_plotdata.csv")
N = len(df)

theta = df["theta"].values
D_obs = df["D_obs"].values
D_pred = df["D_pred"].values
n_obs = df["n_obs"].values
n_pred = df["n_BE_pred"].values
gbar_log = df["gbar_log"].values

resid = np.log10(D_obs) - np.log10(D_pred)
resid_std = np.std(resid)
resid_mean = np.mean(resid)

ln2 = np.log(2.0)
theta_grid = np.linspace(0.08, theta.max(), 600)
D_curve = 1.0 + 1.0 / (np.exp(theta_grid) - 1.0)
n_curve = 1.0 / (np.exp(theta_grid) - 1.0)

# Binned medians for panel (c)
bin_edges = np.logspace(np.log10(theta.min()), np.log10(theta.max()), 12)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
bin_n_med = np.full(len(bin_centers), np.nan)
bin_n_lo = np.full(len(bin_centers), np.nan)
bin_n_hi = np.full(len(bin_centers), np.nan)
for i in range(len(bin_centers)):
    mask = (theta >= bin_edges[i]) & (theta < bin_edges[i+1])
    if mask.sum() >= 5:
        vals = n_obs[mask]
        bin_n_med[i] = np.median(vals)
        bin_n_lo[i] = np.percentile(vals, 16)
        bin_n_hi[i] = np.percentile(vals, 84)

fig = plt.figure(figsize=(7.2, 7.5))
gs = GridSpec(3, 1, height_ratios=[1.35, 1.0, 1.35], hspace=0.75)

# Panel (a): D vs theta with BE curve
ax1 = fig.add_subplot(gs[0])
ax1.scatter(theta, D_obs, s=3, alpha=0.25, color="#1f77b4",
            rasterized=True, label=f"SPARC RAR ($N={N}$)")
ax1.plot(theta_grid, D_curve, color="black", lw=1.8,
         label=r"$D = 1 + n_{\mathrm{BE}}(\theta) = 1/(1 - e^{-\theta})$")
ax1.axvline(ln2, color="crimson", lw=1.0, ls="--", alpha=0.8)
ax1.axhline(2.0, color="crimson", lw=1.0, ls="--", alpha=0.8)
ax1.scatter([ln2], [2.0], s=40, color="crimson", zorder=5,
            edgecolor="black", linewidth=0.6,
            label=r"Crossover $(\ln 2,\,2)$")
ax1.set_xlabel(r"$\theta = \sqrt{g_{\mathrm{bar}}/a_0}$")
ax1.set_ylabel(r"$D = g_{\mathrm{obs}}/g_{\mathrm{bar}}$")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim(0.08, 8)
ax1.set_ylim(0.7, 60)
ax1.legend(loc="upper right", frameon=True, framealpha=0.92)
ax1.set_title(r"(a) Mass discrepancy across the RAR")

# Panel (b): residuals
ax2 = fig.add_subplot(gs[1])
ax2.scatter(gbar_log, resid, s=3, alpha=0.3, color="#1f77b4", rasterized=True)
ax2.axhline(0, color="black", lw=1.0)
ax2.axhline(resid_std, color="crimson", lw=0.9, ls="--", alpha=0.9,
            label=fr"$\pm\sigma = \pm{resid_std:.3f}$ dex")
ax2.axhline(-resid_std, color="crimson", lw=0.9, ls="--", alpha=0.9)
ax2.set_xlabel(r"$\log_{10}(g_{\mathrm{bar}}/[\mathrm{m\,s^{-2}}])$")
ax2.set_ylabel(r"$\log_{10} D_{\mathrm{obs}} - \log_{10} D_{\mathrm{pred}}$")
ax2.set_xlim(gbar_log.min() - 0.1, gbar_log.max() + 0.1)
ax2.set_ylim(-0.8, 0.8)
ax2.legend(loc="upper right", frameon=True, framealpha=0.92)
ax2.set_title(fr"(b) Residuals: mean $={resid_mean:+.3f}$ dex, $\sigma = 0.133$ dex "
              r"$\approx$ intrinsic RAR scatter")

# Panel (c): n_BE predicted vs observed (binned)
ax3 = fig.add_subplot(gs[2])
ax3.scatter(theta, n_obs, s=3, alpha=0.20, color="#1f77b4", rasterized=True)
valid = ~np.isnan(bin_n_med)
ax3.errorbar(bin_centers[valid], bin_n_med[valid],
             yerr=[bin_n_med[valid] - bin_n_lo[valid],
                   bin_n_hi[valid] - bin_n_med[valid]],
             fmt="s", color="#ff7f0e", markersize=5,
             capsize=2.5, elinewidth=0.9, label=r"Binned median $\pm$ 16/84 %")
ax3.plot(theta_grid, n_curve, color="black", lw=1.8,
         label=r"$n_{\mathrm{BE}}(\theta) = 1/(e^\theta - 1)$")
ax3.axvline(ln2, color="crimson", lw=1.0, ls="--", alpha=0.8)
ax3.axhline(1.0, color="crimson", lw=1.0, ls="--", alpha=0.8)
ax3.scatter([ln2], [1.0], s=40, color="crimson", zorder=5,
            edgecolor="black", linewidth=0.6,
            label=r"One bit: $(\ln 2,\,1)$")
ax3.set_xlabel(r"$\theta = \sqrt{g_{\mathrm{bar}}/a_0}$")
ax3.set_ylabel(r"$n = D - 1$")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(0.08, 8)
ax3.set_ylim(0.01, 60)
ax3.legend(loc="upper right", frameon=True, framealpha=0.92, ncol=1)
ax3.set_title(r"(c) Occupation number across the full $\theta$ range")

fig.savefig("mond_be_rar_figure.pdf", dpi=300, bbox_inches="tight")
fig.savefig("mond_be_rar_figure.png", dpi=200, bbox_inches="tight")
print("Figure saved. Summary:")
print(f"  N = {N}")
print(f"  residual mean = {resid_mean:+.5f} dex")
print(f"  residual std  = {resid_std:.5f} dex")
print(f"  (ln 2)^2 = {ln2**2:.5f}")
