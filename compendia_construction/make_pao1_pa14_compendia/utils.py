"""
Author: Alexandra Lee
Date Created: 17 April 2020

Scripts to process expression data for pilot experiment
"""
import pandas as pd
import random
import gzip
from Bio import SeqIO


def get_pao1_pa14_gene_map(gene_annotation_filename, reference_genotype):
    """
    Returns file mapping PAO1 gene ids to PA14 ids and label which genes are core and
    the number of genes that are mapped

    Arguments
    ----------
    gene_annotation_file: str
        File containing mapping between PAO1 and PA14 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PAO1_ID,Name,Product.Name,PA14_ID

    reference_genotype: str
        Either 'pao1' or 'pa14'

    """
    # Check that genotype is set correctly
    assert reference_genotype.lower() in [
        "pao1",
        "pa14",
    ], "Reference genotype string needs to be either pao1 or pa14"

    # Read gene annotation
    gene_mapping = pd.read_csv(gene_annotation_filename, sep=",", header=0)

    reference_id = "PAO1_ID" if reference_genotype.lower() == "pao1" else "PA14_ID"
    non_reference_id = "PA14_ID" if reference_genotype.lower() == "pao1" else "PAO1_ID"

    # Accessory are genes that don't have a homologous gene
    unmapped_genes = gene_mapping[gene_mapping[non_reference_id].isna()]
    acc_genes = list(unmapped_genes[reference_id])

    # Add label for core genes
    gene_mapping.loc[~gene_mapping[reference_id].isin(acc_genes), "annotation"] = "core"

    # Add column with number of genes mapped in case user
    # would like to only consider genes with 1-1 mapping
    gene_mapping["num_mapped_genes"] = (
        gene_mapping[non_reference_id].str.split(", ").str.len()
    )

    gene_mapping.set_index(reference_id, inplace=True)

    return gene_mapping


def get_core_genes(pao1_ref_mapping_df, pa14_ref_mapping_df, is_mapping_1_1):
    """
    Returns list of core genes using PAO1 ids and PA14 ids

    Arguments
    ----------
    pao1_ref_mapping_df: df
        Dataframe generated from get_pao1_pa14_gene_map() to give mapping
        from PAO1 ids to PA14 ids

        Columns: PAO1_ID, Name, Product.Name, PA14_ID, annotation, num_mapped_genes

    pa14_ref_mapping_df: df
        Dataframe generated from get_pao1_pa14_gene_map() to give mapping
        from PA14 ids to PAO1 ids

        Columns: PA14_ID, Name, Product.Name, PAO1_ID, annotation, num_mapped_genes

    is_mapping_1_1: bool
        True if only want to return core genes that have a 1-1 mapping between PAO1 and PA14
    """

    if is_mapping_1_1:
        # Include only core genes that have a 1-1 gene mapping
        pao1_ref_core_df = pao1_ref_mapping_df[
            (pao1_ref_mapping_df["annotation"] == "core")
            & (pao1_ref_mapping_df["num_mapped_genes"] == 1.0)
        ]

        pa14_ref_core_df = pa14_ref_mapping_df[
            (pa14_ref_mapping_df["annotation"] == "core")
            & (pa14_ref_mapping_df["num_mapped_genes"] == 1.0)
        ]
    else:
        pao1_ref_core_df = pao1_ref_mapping_df[
            (pao1_ref_mapping_df["annotation"] == "core")
        ]
        pa14_ref_core_df = pa14_ref_mapping_df[
            (pa14_ref_mapping_df["annotation"] == "core")
        ]

    # Get list of pao1 core genes
    pao1_core_genes = pd.DataFrame(
        list(pao1_ref_core_df.index) + list(pa14_ref_core_df["PAO1_ID"].values),
        columns=["gene id"],
    )
    # Reshape to get single value per row
    pao1_core_genes = pd.DataFrame(
        pao1_core_genes["gene id"].str.split(", ").sum(), columns=["gene id"]
    )
    # Remove duplicates that might exist after taking the union
    pao1_core_genes.drop_duplicates(keep="first", inplace=True)
    pao1_core_genes = list(pao1_core_genes["gene id"])

    # Get list of pa14 core genes
    pa14_core_genes = pd.DataFrame(
        list(pa14_ref_core_df.index) + list(pao1_ref_core_df["PA14_ID"].values),
        columns=["gene id"],
    )
    # Reshape to get single value per row
    pa14_core_genes = pd.DataFrame(
        pa14_core_genes["gene id"].str.split(", ").sum(), columns=["gene id"]
    )
    # Remove duplicates that might exist after taking the union
    pa14_core_genes.drop_duplicates(keep="first", inplace=True)
    pa14_core_genes = list(pa14_core_genes["gene id"])

    return pao1_core_genes, pa14_core_genes


def get_my_core_acc_genes(pao1_annotation_filename, pa14_annotation_filename, my_expression_pao1_df, my_expression_pa14_df):
    """
    Returns lists of core and accessory genes in your dataset
    using both PAO1 and PA14 reference genotypes

    pao1_annotation_file: str
        File containing mapping between PAO1 and PA14 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PAO1_ID,Name,Product.Name,PA14_ID

    pa14_annotation_file: str
        File containing mapping between PA14 and PAO1 gene ids downloaded from
        BACTOME annotation tool: https://pseudomonas-annotator.shinyapps.io/pa_annotator/

        Columns: PA14_ID,Name,Product.Name,PAO1_ID

    my_expression_pao1_df: dataframe
        Gene expression dataframe aligned using PAO1 reference

    my_expression_pa14_df: dataframe
        Gene expression dataframe aligned using PA14 reference
    """
    # Get mapping between PAO1 and PA14 genes using PAO1 and PA14 references
    gene_mapping_pao1 = get_pao1_pa14_gene_map(pao1_annotation_filename, "pao1")
    gene_mapping_pa14 = get_pao1_pa14_gene_map(pa14_annotation_filename, "pa14")

    # Get core genes: genes that have a homolog between PAO1 and PA14
    core_pao1_genes, core_pa14_genes = get_core_genes(
        gene_mapping_pao1,
        gene_mapping_pa14,
        False
    )

    print(f"Number of PAO1 core genes: {len(core_pao1_genes)}")
    print(f"Number of PA14 core genes: {len(core_pa14_genes)}")

    # Select only core genes that are included in my dataset
    pao1_ref_genes = my_expression_pao1_df.columns
    my_core_pao1_genes = list(set(core_pao1_genes).intersection(pao1_ref_genes))

    print(f"Number of PAO1 core genes in my dataset: {len(my_core_pao1_genes)}")

    # Select only core genes that are included in my dataset
    pa14_ref_genes = my_expression_pa14_df.columns
    my_core_pa14_genes = list(set(core_pa14_genes).intersection(pa14_ref_genes))

    print(f"Number of PA14 core genes in my dataset: {len(my_core_pa14_genes)}")

    # Get PAO1-specific genes
    pao1_acc = list(set(pao1_ref_genes) - set(my_core_pao1_genes))
    print(f"Number of PAO1-specific genes: {len(pao1_acc)}")

    # Get PA14-specific genes
    pa14_acc = list(set(pa14_ref_genes) - set(my_core_pa14_genes))
    print(f"Number of PA14-specific genes: {len(pa14_acc)}")

    return {"core_pao1": my_core_pao1_genes,
            "core_pa14": my_core_pa14_genes,
            "acc_pao1": pao1_acc,
            "acc_pa14": pa14_acc
            }
