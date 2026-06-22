# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def _concat_df(files, cols, label):
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        sub = df[cols].copy()
        sub["y"] = label
        frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=cols + ["y"])
    return pd.concat(frames, ignore_index=True)


cols = ["Danc", "Dplus", "df", "fd", "DD", "U", "Q"]
stats = {
    "Danc": "$D_{anc}$",
    "Dplus": "$D^+$",
    "df": "$d_f$",
    "fd": "$f_d$",
    "DD": "$D_D$",
    "U": "$U50$",
    "Q": "$Q95$",
}

neu_all = list(snakemake.input.neu)
sel_all = list(snakemake.input.sel)

np.random.seed(int(snakemake.params.seed))
R = int(snakemake.params.replicates)
seed_list = np.random.randint(1, 2**31, R)

fpr_grid = np.linspace(0, 1, 200)
mean_curves = {col: [] for col in cols}
auc_logs = {col: [] for col in cols}

for r in range(R):
    rng = np.random.default_rng(seed_list[r])
    neu_files = rng.choice(neu_all, 500, replace=False)
    sel_files = rng.choice(sel_all, 500, replace=False)

    df_neu = _concat_df(neu_files, cols=cols, label=0)
    df_sel = _concat_df(sel_files, cols=cols, label=1)
    df = pd.concat([df_neu, df_sel], ignore_index=True)
    y_true = df["y"].to_numpy().copy()

    for col in cols:
        y_score = df[col].to_numpy().copy()
        mask = np.isfinite(y_score)
        if mask.sum() == 0:
            continue
        fpr, tpr, _ = roc_curve(y_true[mask], y_score[mask])
        roc_auc = auc(fpr, tpr)
        auc_logs[col].append(roc_auc)
        tpr_interp = np.interp(fpr_grid, fpr, tpr)
        mean_curves[col].append(tpr_interp)

plt.figure(figsize=(5, 5), dpi=300)

for col in cols:
    if not auc_logs[col]:
        continue
    arr_auc = np.array(auc_logs[col])
    mean_auc = arr_auc.mean()
    mean_fpr = np.mean(mean_curves[col], axis=0)
    plt.plot(fpr_grid, mean_fpr, lw=2, label=f"{stats[col]} (mean AUC = {mean_auc:.2f})")

plt.plot([0, 1], [0, 1], "--", lw=1)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"$s$ = {snakemake.wildcards.s}")
plt.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
