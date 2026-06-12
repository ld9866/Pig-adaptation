#!/usr/bin/env python3
import os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

TPM_FILE    = "TPM_matrix.csv"
GROUP_FILE  = "sample_group.tsv"                       
TARGET_FILE = "target_genes.txt"                       
GTF_FILE    = "00.index/Sus_scrofa.Sscrofa11.1.111.chr.gtf"  
OUT_PDF     = "gene_boxplots.pdf"

tpm    = pd.read_csv(TPM_FILE, index_col=0)
groups = pd.read_csv(GROUP_FILE, sep="\t").set_index("Sample")["Group"]
targets = [t.strip() for t in open(TARGET_FILE) if t.strip()]

name2id = {}
if os.path.exists(GTF_FILE):
    with open(GTF_FILE) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            gnm = re.search(r'gene_name "([^"]+)"', f[8])
            if gid and gnm:
                name2id[gnm.group(1)] = gid.group(1)

resolved = []
for t in targets:
    if t in tpm.index:
        resolved.append((t, t))
    elif name2id.get(t) in tpm.index:
        resolved.append((t, name2id[t]))
    else:
        print("skip:", t)

sub = tpm.loc[[rid for _, rid in resolved]].T
sub.columns = [lab for lab, _ in resolved]
sub["Group"] = groups.reindex(sub.index).values

sns.set_theme(style="whitegrid")
labels = [lab for lab, _ in resolved]
order  = sorted(sub["Group"].dropna().unique())
pal    = dict(zip(order, sns.color_palette("Set2", len(order))))

ncol = min(4, max(1, len(labels)))
nrow = (len(labels) + ncol - 1) // ncol
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow), squeeze=False)
axes = axes.flatten()

for i, g in enumerate(labels):
    d = sub[["Group", g]].dropna()
    sns.boxplot(data=d, x="Group", y=g, hue="Group", order=order,
                palette=pal, dodge=False, width=0.6, ax=axes[i])
    if axes[i].get_legend():
        axes[i].get_legend().remove()
    sns.stripplot(data=d, x="Group", y=g, order=order,
                  color=".3", size=4, alpha=0.7, ax=axes[i])
    axes[i].set_title(g, fontweight="bold")
    axes[i].set_xlabel("")
    axes[i].set_ylabel("TPM")
for j in range(len(labels), len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(OUT_PDF, dpi=300)
print("save", OUT_PDF, "(", len(labels), "gene )")
