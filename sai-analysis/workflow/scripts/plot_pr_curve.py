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
from sklearn.metrics import precision_recall_curve, average_precision_score


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

wo_all = list(snakemake.input.wo)
neu_all = list(snakemake.input.neu)
sel_all = list(snakemake.input.sel)

recall_grid = [0, 0.5, 1]
mean_curves = {col: [] for col in cols}
ap_logs = {col: [] for col in cols}

np.random.seed(int(snakemake.params.seed))
R = int(snakemake.params.replicates)
seed_list = np.random.randint(1, 2**31, R)

for r in range(R):
    rng = np.random.default_rng(seed_list[r])
    wo_files = rng.choice(wo_all, 980, replace=False)
    neu_files = rng.choice(neu_all, 18, replace=False)
    sel_files = rng.choice(sel_all, 2, replace=False)

    df_wo = _concat_df(wo_files, cols=cols, label=0)
    df_neu = _concat_df(neu_files, cols=cols, label=0)
    df_sel = _concat_df(sel_files, cols=cols, label=1)
    df = pd.concat([df_wo, df_neu, df_sel], ignore_index=True)
    y_true = df["y"].to_numpy().copy()

    for col in cols:
        y_score = df[col].to_numpy().copy()
        mask = np.isfinite(y_score)
        if mask.sum() == 0 or y_true[mask].sum() != 2: # ensure two positive instances' scores are not nan
            continue
        precision, recall, _ = precision_recall_curve(y_true[mask], y_score[mask])
        ap = average_precision_score(y_true[mask], y_score[mask])
        ap_logs[col].append(ap)

        uniq_recall, idx = np.unique(recall, return_inverse=True)
        # uniq_recall: [0, 0.5, 1]
        max_precision = np.full_like(uniq_recall, -np.inf, dtype=float)
        np.maximum.at(max_precision, idx, precision)

        if np.unique((y_score[mask][-2:])).size == 1:
            # uniq_recall: [0, 1]
            max_precision = np.insert(max_precision.astype(float), 1, max_precision[1])

        mean_curves[col].append(max_precision)

plt.figure(figsize=(5, 5), dpi=300)

pos_rate = 0.002
plt.hlines(pos_rate, 0, 1, linestyles="--", linewidth=1, label=f"baseline = {pos_rate:.3f}")

for col in cols:
    if not ap_logs[col]:
        continue
    arr_ap = np.array(ap_logs[col])
    mean_ap = arr_ap.mean()
    mean_prec = np.mean(mean_curves[col], axis=0)
    # https://scikit-learn.org/0.19/auto_examples/model_selection/plot_precision_recall.html#plot-the-precision-recall-curve
    plt.step(recall_grid[::-1], mean_prec[::-1], where="post", lw=2, label=f"{stats[col]} (mean AP = {mean_ap:.2f})")

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title(f"$s$ = {snakemake.wildcards.s}")
plt.legend(loc="lower left", fontsize=8)
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
