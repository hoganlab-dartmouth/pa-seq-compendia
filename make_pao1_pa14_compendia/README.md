# Make PAO1 and PA14 compendia

These notebooks aim to partition the heuristic filtered gene expression data into a PAO1 and PA14 compendia where the PAO1 compendium contains only PAO1-like samples aligned to the PAO1 reference and the PA14 compendium contains only PA14-like samples aligned to the PA14 reference. 
To determine which samples are PAO1 versus PA14 we will use the median expression of accessory genes to determine if a sample is PAO1 or PA14.
In our exploratory analysis (`cluster_by_accessory_gene.ipynb`) we found that samples labeled as PAO1 based on SRA annotations had high PAO1 accessory gene expression.
Whereas samples labeled as PA14 by SRA had high PA14 accessory gene expression.

See plot below where the median expression of PAO1-only genes (PAO1 accessory genes) on the x-axis and the median expression of PA14-only genes (PA14 accessory genes) on the y-axis.
Each point is a sample.
![all_samples](https://github.com/greenelab/core-accessory-interactome/blob/master/explore_data/Expression_accessory_genes_all_samples.svg)

## Directory Structure
| Folder | Description |
| --- | --- |
| [0_decide_threshold.ipynb](0_decide_threshold.ipynb) | This analysis notebook examines the distribution of accessory gene expression to determine the threshold that will be used to partition the data.|
| [1_create_compendia.ipynb](1_create_compendia.ipynb) | This analysis notebook uses the thresholds determined in the previous notebook to partition the gene expression data into PAO1 and PA14 compendia.|
| [2_validate_compendia.ipynb](2_validate_compendia.ipynb) | This analysis notebook visualizes the compendia to make sure the data looks as expected. These visualizations can be found in figure 3.|

## Dependencies

We use Python and Jupyter notebooks to load and process the data. 
These analyses use a different environment from the rest of the repository.
The environment dependencies are specificed in the `environment.yml` file in this directory and can be installed by running the following terminal commands from this directory:
```
conda env create -f environment.yml
conda activate make_compendia
```