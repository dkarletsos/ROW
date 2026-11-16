import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
INK, INK2, GRID, ROW = "#1B1F24", "#4E5560", "#E3E6EA", "#5E4FC9"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"axes.edgecolor":INK2,"xtick.color":INK2,"ytick.color":INK2})
R={"ROW":(1.29,1.40,1.64),"AIPW-T1/99":(1.31,1.40,1.64),"AIPW-T5/95":(1.71,1.38,1.61),"AEB-main":(1.32,1.44,1.78),"SU":(2.18,1.38,1.60),
"G-computation":(2.36,1.37,1.60),"TMLE":(2.03,1.98,2.36),"EB-main":(2.14,1.23,2.70),"AIPW":(2.39,1.80,2.72),"AOW":(3.08,1.81,1.85),"OW":(3.13,1.87,2.29)}
order=sorted(R,key=lambda m:(max(R[m]),m!="ROW"))
fig,ax=plt.subplots(figsize=(7.2,5.0))
mk=[("o","#8C93A0","Main design (96)"),("s","#C9A227","Held-out design (48)"),("^","#2E8B8B","Stress-test design (128)")]
for i,m in enumerate(order):
    v=R[m]; col=ROW if m=="ROW" else INK2
    ax.plot([min(v),max(v)],[i,i],color=GRID,lw=6,solid_capstyle="round",zorder=1)
    for j,(s,c,l) in enumerate(mk): ax.scatter(v[j],i,marker=s,s=46,color=c,edgecolor="white",linewidth=0.8,zorder=3,label=l if i==0 else None)
    ax.text(max(v)+0.05,i,f"{max(v):.2f}",va="center",fontsize=8.5,color=col,fontweight="bold" if m=="ROW" else "normal")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
for t in ax.get_yticklabels():
    if t.get_text()=="ROW": t.set_color(ROW); t.set_fontweight("bold")
ax.invert_yaxis(); ax.axvline(1,color=INK2,lw=0.8,ls=":"); ax.set_xlim(0.95,3.35); ax.grid(axis="x",color=GRID); ax.set_axisbelow(True)
ax.set_xlabel("Worst-case regret (RMSE relative to best ATE estimator in each combination)")
ax.legend(frameon=False,fontsize=8,loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=3); ax.set_title("Worst-case regret by design; label = worst over all 272 combinations",fontsize=9.5,color=INK,loc="left")
fig.tight_layout(); fig.savefig("/home/claude/fig2_worst_case.png",dpi=220,facecolor="white")
