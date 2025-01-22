import sys
sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
pd.options.mode.chained_assignment = None

def get_performance(data, geno_state):
    sub_data = data[data['Genotype'] == geno_state]
    sub_data['Precision'] = sub_data['True_positives_length'] / sub_data['Inferred_tracts_length'] * 100
    sub_data['Recall'] = sub_data['True_positives_length'] / sub_data['True_tracts_length'] * 100
    mean_performance = sub_data.groupby("Cutoff").mean(numeric_only=True).dropna()[['Precision', 'Recall']]
    var_performance = sub_data.groupby("Cutoff").var(numeric_only=True).dropna()[['Precision', 'Recall']]
    std_performance = np.sqrt(var_performance)

    return mean_performance, std_performance

data = pd.read_csv(snakemake.input.performance, sep="\t")

phased_data_mean_performance, phased_data_std_performance = get_performance(data, 'phased')
unphased_data_mean_performance, unphased_data_std_performance = get_performance(data, 'unphased')

plt.errorbar(phased_data_mean_performance['Recall'],
             phased_data_mean_performance['Precision'], marker='o', label='Phased data',
             xerr=phased_data_std_performance['Recall'], yerr=phased_data_std_performance['Precision'])
plt.errorbar(unphased_data_mean_performance['Recall'],
             unphased_data_mean_performance['Precision'], marker='o', label='Unphased data',
             xerr=unphased_data_std_performance['Recall'], yerr=unphased_data_std_performance['Precision'])
plt.plot([0,100],[2,2], label='Baseline', linestyle='dashed')
plt.xlim([0,100])
plt.ylim([0,100])
plt.title('PR curve')
plt.xlabel('Recall (%)')
plt.ylabel('Precision (%)')
plt.legend(loc="lower left")

plt.savefig(snakemake.output.pdf, bbox_inches='tight')
