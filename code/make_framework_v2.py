import sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "figures"
INK, INK2 = "#0b0b0b", "#52514e"
Q, A_ATE, A_ATO, A_WARN = "#e8eef6", "#d6e6f8", "#fbe1d6", "#f3f3f1"

fig, ax = plt.subplots(figsize=(10, 5.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")

def box(x, y, w, h, text, fc, bold_first=True, fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                fc=fc, ec=INK2, lw=0.8))
    lines = text.split("\n")
    ax.text(x + w / 2, y + h - 1.6, lines[0], ha="center", va="top", fontsize=fs + 0.6,
            fontweight="bold" if bold_first else "normal", color=INK)
    if len(lines) > 1:
        ax.text(x + w / 2, y + h - 4.4, "\n".join(lines[1:]), ha="center", va="top",
                fontsize=fs, color=INK, linespacing=1.35)

def arrow(x1, y1, x2, y2, label=None, lx=0, ly=0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11, color=INK2, lw=1))
    if label: ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, fontsize=8, color=INK2, ha="center",
                      style="italic")

box(26, 47, 48, 8, "Step 1  Define the target population (ICH E9(R1))\nWho should the comparison describe?", Q)
arrow(40, 47, 20, 40.5, "all eligible patients", -7, 0)
arrow(60, 47, 80, 40.5, "patients in clinical equipoise", 9, 0)

box(1, 30, 38, 10, "Population-level question  →  ATE\ne.g. 'what if every eligible patient\nreceived each option?'", A_ATE)
box(61, 30, 38, 10, "Equipoise question  →  ATO\ne.g. 'among patients who could realistically\nreceive any option, which is better?'", A_ATO)

box(30, 17.2, 40, 10.3, "Step 2  Diagnose overlap\nGPS densities by arm · share with min GPS < 0.05\nESS under untrimmed IPTW · max |SMD|", Q)
arrow(20, 30, 36, 27.5); arrow(80, 30, 64, 27.5)

box(1, 1, 30, 14, "Good overlap\nESS ≥ ~80%, few extreme GPS\n\nATE via IPTW;\ntrimming adds little.", A_WARN)
box(35, 1, 30, 14, "Poor overlap, ATE required\nDR-ROW + flexible outcome model.\nFlat λ path → small λ;\ndrifting → improve the model.\nNo trusted model: ROW, λ set in SAP.", A_WARN)
box(69, 1, 30, 14, "Poor overlap, ATO acceptable\n\nOverlap weights: lowest RMSE,\nnominal coverage, best balance.\nLabel it ATO; describe population.", A_WARN)
arrow(40, 16.6, 16, 15.4); arrow(50, 16.6, 50, 15.4); arrow(60, 16.6, 84, 15.4)

fig.savefig(f"{OUT}/fig_decision_framework_v2.png", dpi=220, bbox_inches="tight", facecolor="white")
print("framework written")
