import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RD = sys.argv[1] if len(sys.argv) > 1 else "results_row"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ROWC, T99, T95 = "#6250d6", "#1baf7a", "#eda100"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
a = pd.read_csv(f"{RD}/row_performance_by_scenario.csv")
a = a[(a.effect == "het") & (a.n == 2000)]
path = ["OW", "ROW-c0.5", "ROW-c2", "ROW-c8", "ROW-c32", "IPTW"]
plab = ["OW\n(λ=0)", "0.5", "2", "8", "32", "IPTW\n(λ=∞)"]
panels = [(3, "poor", "none"), (3, "good", "severe"), (4, "poor", "none"), (4, "poor", "severe")]
rows = [("rmse_ATE", "RMSE vs true ATE", None), ("cov_ATE", "Coverage of ATE (%)", 95), ("cov", "Coverage of own target (%)", 95)]
fig, ax = plt.subplots(3, 4, figsize=(11, 7.2), sharex=True)
for j, (K, ov, pos) in enumerate(panels):
    d = a[(a.K == K) & (a.overlap == ov) & (a.positivity == pos)].set_index("method")
    for i, (col, lab, ref) in enumerate(rows):
        x = ax[i, j]; x.grid(axis="y", color=GRID); x.set_axisbelow(True)
        sc = 100 if col.startswith("cov") else 1
        y = [d.loc[m, col] * sc for m in path]
        x.plot(range(6), y, color=ROWC, marker="o", lw=2, ms=6, label="ROW path (OW → IPTW)")
        if col != "cov":
            x.axhline(d.loc["IPTW-T1/99", col] * sc, color=T99, lw=1.8, ls="--", label="IPTW trim 1/99")
            x.axhline(d.loc["IPTW-T5/95", col] * sc, color=T95, lw=1.8, ls=":", label="IPTW trim 5/95")
        if ref: x.axhline(ref, color=INK2, lw=0.8)
        if col.startswith("cov"): x.set_ylim(0, 100)
        if i == 0: x.set_title(f"{K} arms, {ov} overlap,\npositivity {'violated' if pos == 'severe' else 'respected'}", fontsize=9)
        if j == 0: x.set_ylabel(lab)
        x.set_xticks(range(6)); x.set_xticklabels(plab, fontsize=8)
    for i in range(3):
        if rows[i][0] == "rmse_ATE":
            hi = max(xx.get_ylim()[1] for xx in ax[i]); [xx.set_ylim(0, hi) for xx in ax[i]]
fig.supxlabel("λ / K²  (regularisation strength)", fontsize=9)
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(f"{OUT}/fig_row_path.png", dpi=200); print("ok")
