# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python [conda env:core_acc]
#     language: python
#     name: conda-env-core_acc-py
# ---

# # Decide threshold
#
# The goal of this notebook is to define a threshold to determine if a sample if PAO1 or not (likewise, if a sample is PA14 or not). We used known labels from SRA to do this. Specifically, we examined the distribution of PAO1 samples (grey) vs non-PAO1 samples (blue). We define the threshold to be one that separated between the two distributions. We use this threshold in [1_create_compendia.ipynb](1_create_compendia.ipynb) to partition gene expression data into PAO1 and PA14 compendia.

# %load_ext autoreload
# %autoreload 2
# %matplotlib inline
import os
import pandas as pd
import seaborn as sns
from textwrap import fill
import matplotlib.pyplot as plt
import utils

# +
# Raw (TPM) expression data files
pao1_expression_filename = "../qc_filtering/qc-out/pao1_aligned_compendium_p2_filtered_counts_norm.csv" 
pa14_expression_filename = "../qc_filtering/qc-out/pa14_aligned_compendium_p2_filtered_counts_norm.csv"

# File containing table to map sample id to strain name
sample_to_strain_filename = "Run_Table_Strain_Bool_GD.csv"
# -

# Load expression data
# Matrices will be sample x gene after taking the transpose
pao1_expression = pd.read_csv(pao1_expression_filename, index_col=0, header=0).T
pa14_expression = pd.read_csv(pa14_expression_filename, index_col=0, header=0).T

# Load metadata
# Set index to experiment id, which is what we will use to map to expression data
sample_to_strain_table_full = pd.read_csv(sample_to_strain_filename, index_col=2)

# ## Format expression data
#
# Format index to only include experiment id. This will be used to map to expression data and SRA labels later

# +
# Format expression data indices so that values can be mapped to `sample_to_strain_table`
pao1_index_processed = pao1_expression.index.str.split(".").str[0]
pa14_index_processed = pa14_expression.index.str.split(".").str[0]

print(
    f"No. of samples processed using PAO1 reference after filtering: {pao1_expression.shape}"
)
print(
    f"No. of samples processed using PA14 reference after filtering: {pa14_expression.shape}"
)
pao1_expression.index = pao1_index_processed
pa14_expression.index = pa14_index_processed
# -

pao1_expression.head()

pa14_expression.head()

# ## Format SRA annotations

# +
# Since experiments have multiple runs there are duplicated experiment ids in the index
# We will need to remove these so that the count calculations are accurate
sample_to_strain_table_full_processed = sample_to_strain_table_full[
    ~sample_to_strain_table_full.index.duplicated(keep="first")
]

assert (
    len(sample_to_strain_table_full.index.unique())
    == sample_to_strain_table_full_processed.shape[0]
)

# +
# Aggregate boolean labels into a single strain label
aggregated_label = []
for exp_id in list(sample_to_strain_table_full_processed.index):
    if sample_to_strain_table_full_processed.loc[exp_id, "PAO1"].all() == True:
        aggregated_label.append("PAO1")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PA14"].all() == True:
        aggregated_label.append("PA14")
    elif sample_to_strain_table_full_processed.loc[exp_id, "PAK"].all() == True:
        aggregated_label.append("PAK")
    elif (
        sample_to_strain_table_full_processed.loc[exp_id, "ClinicalIsolate"].all()
        == True
    ):
        aggregated_label.append("Clinical Isolate")
    else:
        aggregated_label.append("NA")

sample_to_strain_table_full_processed["Strain type"] = aggregated_label

sample_to_strain_table = sample_to_strain_table_full_processed["Strain type"].to_frame()

sample_to_strain_table.head()
# -

# ## Label samples with SRA annotations

# Label samples with SRA annotations
# pao1_expression_label = pao1_expression_binned.join(
#    sample_to_strain_table, how='left')
pao1_expression_label = pao1_expression.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
pa14_expression_label = pa14_expression.merge(
    sample_to_strain_table, left_index=True, right_index=True
)
print(pao1_expression_label.shape)
pao1_expression_label.head()

print(pa14_expression_label.shape)
pa14_expression_label.head()

# ## Get accessory gene expression

# +
pao1_annot_filename = "PAO1_ID_2_PA14_ID_PAO1ref.csv" 
pa14_annot_filename = "PA14_ID_2_PAO1_ID_PA14ref.csv"

core_acc_dict = utils.get_my_core_acc_genes(
    pao1_annot_filename, pa14_annot_filename, pao1_expression, pa14_expression
)
# -

pao1_acc = core_acc_dict["acc_pao1"]
pa14_acc = core_acc_dict["acc_pa14"]

# +
# Create accessory df
# accessory gene ids | median accessory expression | strain label

# PAO1
pao1_acc_expression = pao1_expression_label[pao1_acc]
pao1_acc_expression["median_acc_expression"] = pao1_acc_expression.median(axis=1)

# PA14
pa14_acc_expression = pa14_expression_label[pa14_acc]
pa14_acc_expression["median_acc_expression"] = pa14_acc_expression.median(axis=1)
# -

# Add back labels
pao1_acc_expression["Strain type"] = pao1_expression_label["Strain type"]
pa14_acc_expression["Strain type"] = pa14_expression_label["Strain type"]

pao1_acc_expression.head()

pa14_acc_expression.head()

# +
# Merge PAO1 and PA14 accessory dataframes
pao1_pa14_acc_expression = pao1_acc_expression.merge(
    pa14_acc_expression,
    left_index=True,
    right_index=True,
    suffixes=["_pao1", "_pa14"],
)

pao1_pa14_acc_expression.head()
# -

# ## Plot distribution

# +
# Get PAO1 samples that are labeled PAO1 and non-PAO1
pao1_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pao1"] == "PAO1",
    "median_acc_expression_pao1",
]

non_pao1_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pao1"] != "PAO1",
    "median_acc_expression_pao1",
]

# Format data for plotting
pao1_acc_expression = pao1_pa14_acc_expression.copy()
pao1_acc_expression.loc[pao1_sra.index, "label"] = "PAO1"
pao1_acc_expression.loc[non_pao1_sra.index, "label"] = "not PAO1"
# -

pao1_threshold = 25

f = sns.displot(
    data=pao1_acc_expression,
    x="median_acc_expression_pao1",
    hue="label",
    palette={
        "not PAO1": "grey",
        "PAO1": "#C6A9B5"
    },
    aspect=1.5,
    edgecolor=None,
    legend=False
           )
plt.axvline(pao1_threshold,0,100, color="red")
plt.title("PAO1 median accessory expression", family="sans-serif", size=16)
plt.xlabel("median accessory expression", family="sans-serif", size=14)
plt.ylabel("count", family="sans-serif", size=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
f.savefig("pao1_threshold_dist.svg", format="svg", dpi=300)

# +
# Get PA14 samples that are labeled PA14 and non-PA14
pa14_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pa14"] == "PA14",
    "median_acc_expression_pa14",
]

non_pa14_sra = pao1_pa14_acc_expression.loc[
    pao1_pa14_acc_expression["Strain type_pa14"] != "PA14",
    "median_acc_expression_pa14",
]

# Format data for plotting
pa14_acc_expression = pao1_pa14_acc_expression.copy()
pa14_acc_expression.loc[pa14_sra.index, "label"] = "PA14"
pa14_acc_expression.loc[non_pa14_sra.index, "label"] = "not PA14"
# -

pa14_threshold = 25

g = sns.displot(
    data=pa14_acc_expression,
    x="median_acc_expression_pa14",
    hue="label",
    palette={
        "not PA14": "grey",
        "PA14": "#C6A9B5"
    },
    aspect=1.5,
    edgecolor=None,
    legend=False
           )
plt.axvline(pa14_threshold,0,100, color="red")
plt.title("PA14 median accessory expression", family="sans-serif", size=16)
plt.xlabel("median accessory expression", family="sans-serif", size=14)
plt.ylabel("count", family="sans-serif", size=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
g.savefig("pa14_threshold_dist.svg", format="svg", dpi=300)

# **Takeaway:**
# Looks like using a threshold of 25 TPM separates between SRA-annotated PAO1 samples vs non-PAO1 samples. Similarly for PA14. This is the threshold we'll use to bin samples into PAO1 vs PA14 compendia.
