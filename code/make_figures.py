import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures"
import os; os.makedirs(OUT, exist_ok=True)

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"IPTW": "#2a78d6", "IPTW-T1/99": "#1baf7a", "IPTW-T5/95": "#eda100", "OW": "#eb6834"}
MK = {"IPTW": "o", "IPTW-T1/99": "s", "IPTW-T5/95": "^", "OW": "D"}
LAB = {"IPTW": "IPTW (untrimmed)", "IPTW-T1/99": "IPTW trimmed 1/99", "IPTW-T5/95": "IPTW trimmed 5/95", "OW": "Overlap weights"}
METHODS = ["IPTW", "IPTW-T1/99", "IPTW-T5/95", "OW"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})

agg = pd.read_csv(f"{RES}/performance_by_scenario.csv")
con = pd.read_csv(f"{RES}/performance_by_contrast.csv")
tr = pd.read_csv(f"{RES}/true_estimands.csv")
ov = ["good", "moderate", "poor"]
xlab = ["Good", "Moderate", "Poor"]

def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------- Figure 1: performance grid (heterogeneous effects, n = 2000) ----------
def perf_grid(effect, n, fname, metrics):
    fig, axes = plt.subplots(len(metrics), 4, figsize=(10, 2.1 * len(metrics)), sharex=True)
    panels = [(3, "none"), (3, "severe"), (4, "none"), (4, "severe")]
    for j, (K, pos) in enumerate(panels):
        d = agg[(agg.effect == effect) & (agg.n == n) & (agg.K == K) & (agg.positivity == pos)]
        for i, (col, ylab, ref) in enumerate(metrics):
            ax = axes[i, j]; style(ax)
            for m in METHODS:
                dm = d[d.method == m].set_index("overlap").loc[ov]
                y = dm[col] * (100 if col.startswith("coverage") else 1)
                ax.plot(range(3), y, color=COL[m], marker=MK[m], lw=2, ms=6, label=LAB[m])
            if ref is not None: ax.axhline(ref, color=INK2, lw=0.8, ls="--")
            if col.startswith("coverage"): ax.set_ylim(0, 100)
            if i == 0: ax.set_title(f"{K} arms, positivity {'violated' if pos=='severe' else 'respected'}", fontsize=9, color=INK)
            if j == 0: ax.set_ylabel(ylab)
            ax.set_xticks(range(3)); ax.set_xticklabels(xlab)
        for i in range(len(metrics)):
            lo = min(a.get_ylim()[0] for a in axes[i]); hi = max(a.get_ylim()[1] for a in axes[i])
            if metrics[i][0].startswith("coverage"): lo, hi = 0, 100
            for a in axes[i]: a.set_ylim(lo, hi)
    fig.supxlabel("Covariate overlap", fontsize=9, color=INK)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(f"{OUT}/{fname}", dpi=220); plt.close(fig)

perf_grid("het", 2000, "fig_performance_het_n2000.png",
          [("mean_abs_bias", "|Bias| vs own target", 0), ("rmse", "RMSE", None),
           ("coverage", "95% CI coverage (%)", 95), ("ess_pct", "ESS (% of n)", None)])
perf_grid("het", 500, "fig_performance_het_n500.png",
          [("mean_abs_bias", "|Bias| vs own target", 0), ("rmse", "RMSE", None),
           ("coverage", "95% CI coverage (%)", 95), ("ess_pct", "ESS (% of n)", None)])

# ---------- Figure 2: estimand gap (true ATE vs ATO, arm 2 vs arm 1) ----------
fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharey=True)
for j, K in enumerate([3, 4]):
    ax = axes[j]; style(ax)
    for k, (pos, ls) in enumerate([("none", "-"), ("severe", "--")]):
        d = tr[(tr.K == K) & (tr.positivity == pos) & (tr.arm == 2)].set_index("overlap").loc[ov]
        ax.plot(range(3), d.ATO_het, color=COL["OW"], marker="D", lw=2, ls=ls, ms=6,
                mfc=COL["OW"] if pos == "none" else "white", mew=1.5, 
                label=f"ATO, positivity {'violated' if pos=='severe' else 'respected'}")
    ax.axhline(-2.0, color=COL["IPTW"], lw=2, label="ATE (all scenarios)")
    ax.set_xticks(range(3)); ax.set_xticklabels(xlab); ax.set_title(f"{K} arms", fontsize=9)
    if j == 0: ax.set_ylabel("True effect, arm 2 vs arm 1")
    ax.set_xlabel("Covariate overlap")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02), handlelength=3.5)
fig.tight_layout(rect=(0, 0, 1, 0.9)); fig.savefig(f"{OUT}/fig_estimand_gap.png", dpi=220); plt.close(fig)

# ---------- Figure 3: what OW estimates if read as the ATE ----------
fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharey=True)
for j, n in enumerate([500, 2000]):
    ax = axes[j]; style(ax)
    d = agg[(agg.effect == "het") & (agg.n == n) & (agg.K == 4) & (agg.positivity == "severe")]
    for m in METHODS:
        dm = d[d.method == m].set_index("overlap").loc[ov]
        ax.plot(range(3), dm.coverage_vs_ATE * 100, color=COL[m], marker=MK[m], lw=2, label=LAB[m])
    ax.axhline(95, color=INK2, lw=0.8, ls="--")
    ax.set_xticks(range(3)); ax.set_xticklabels(xlab); ax.set_ylim(0, 100)
    ax.set_title(f"n = {n}", fontsize=9)
    if j == 0: ax.set_ylabel("Coverage of the true ATE (%)")
    ax.set_xlabel("Covariate overlap")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
fig.tight_layout(rect=(0, 0, 1, 0.88)); fig.savefig(f"{OUT}/fig_coverage_vs_ATE.png", dpi=220); plt.close(fig)
print("figures written")
