"""Figure: change in RMSE for the ATE, ROW versus overlap weights, n = 2,000 (paired simulation)."""
import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
IN, OUT = (sys.argv[1:3] + ["results_paired", "figures"][len(sys.argv) - 1:])[:2]
INK, INK2, GRID, ROW, ROWL = "#1B1F24", "#4E5560", "#E3E6EA", "#5E4FC9", "#A99FE6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
d = pd.read_csv(f"{IN}/paired_vs_ow.csv"); d = d[(d.method == "ROW") & (d.n == 2000)]
sc = [(K, pos, ov) for K in (3, 4) for pos in ("none", "severe") for ov in ("good", "moderate", "poor")]
lab = [f"{K} arms, {ov}{', PV' if pos == 'severe' else ''}" for K, pos, ov in sc]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.9), sharey=True)
for ax, eff, title, series in [
    (axes[0], "het", "Heterogeneous effects (ATE ≠ ATO)", [("quad", ROW, "Outcome model correct"), ("lin", ROWL, "Outcome model misspecified")]),
    (axes[1], "hom", "Homogeneous effects (ATE = ATO)", [("lin", ROW, "Outcome model correct (linear)"), ("quad", ROWL, "Outcome model correct (flexible)")])]:
    for j, (om, col, nm) in enumerate(series):
        v = [100 * (d[(d.effect == eff) & (d.om == om) & (d.K == K) & (d.positivity == pos) & (d.overlap == ov)].ratio.iloc[0] - 1) for K, pos, ov in sc]
        ax.barh([i + (j - 0.5) * 0.38 for i in range(len(sc))], v, height=0.36, color=col, label=nm)
    ax.axvline(0, color=INK2, lw=0.9); ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    ax.set_xlim(-70, 30); ax.set_title(title, fontsize=10, color=INK); ax.set_xlabel("Change in RMSE versus overlap weights (%)")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
axes[0].set_yticks(range(len(sc))); axes[0].set_yticklabels(lab); axes[0].invert_yaxis()
fig.tight_layout(); fig.savefig(f"{OUT}/fig_row_vs_ow.png", dpi=220, facecolor="white"); print("ok")
