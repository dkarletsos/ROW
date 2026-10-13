import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DD, MD, OUT = (sys.argv[1:4] + ["results_drrow", "results", "figures"][len(sys.argv) - 1:])[:3]
INK2, GRID, LIN, QUAD, T99 = "#52514e", "#e4e3df", "#6250d6", "#e87ba4", "#1baf7a"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
a = pd.read_csv(f"{DD}/drrow_performance_by_scenario.csv")
a = a[(a.n == 2000) & (a.effect == "het") & a.method.str.startswith("DR")]
m = pd.read_csv(f"{MD}/performance_by_scenario.csv")
m = m[(m.n == 2000) & (m.effect == "het") & (m.method == "IPTW-T1/99")]
cs = ["0", "0.5", "2", "8", "32", "Inf"]; lab = ["0\n(OW)", "0.5", "2", "8", "32", "∞\n(AIPW)"]
panels = [(3, "poor", "none"), (3, "poor", "severe"), (4, "poor", "none"), (4, "poor", "severe")]
fig, ax = plt.subplots(2, 4, figsize=(11, 5.4), sharex=True)
for j, (K, ov, pos) in enumerate(panels):
    for i, (col, yl) in enumerate([("rmse", "RMSE vs true ATE"), ("cov", "Coverage of ATE (%)")]):
        x = ax[i, j]; x.grid(axis="y", color=GRID); x.set_axisbelow(True)
        sc = 100 if col == "cov" else 1
        for om, colr, ls, nm in [("lin", LIN, "-", "DR-ROW, linear outcome model (misspecified)"),
                                 ("quad", QUAD, "--", "DR-ROW, flexible outcome model (correct)")]:
            d = a[(a.K == K) & (a.overlap == ov) & (a.positivity == pos) & (a.om == om)].set_index("method")
            x.plot(range(6), [d.loc["DR-ROW-c" + c, col] * sc for c in cs], color=colr, ls=ls, marker="o", lw=2, ms=5, label=nm)
        r = m[(m.K == K) & (m.overlap == ov) & (m.positivity == pos)]
        x.axhline(float(r["rmse" if col == "rmse" else "coverage"].iloc[0]) * sc, color=T99, ls=":", lw=1.8, label="IPTW trim 1/99 (reference)")
        if col == "cov": x.axhline(95, color=INK2, lw=0.8); x.set_ylim(0, 100)
        else: x.set_ylim(0, 0.7)
        if i == 0: x.set_title(f"{K} arms, poor overlap,\npositivity {'violated' if pos == 'severe' else 'respected'}", fontsize=9)
        if j == 0: x.set_ylabel(yl)
        x.set_xticks(range(6)); x.set_xticklabels(lab, fontsize=8)
fig.supxlabel("λ / K²", fontsize=9)
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(f"{OUT}/fig_drrow_path.png", dpi=200); print("ok")
