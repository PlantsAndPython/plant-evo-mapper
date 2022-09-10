# The topological shape of gene expression across the evolution of flowering plants

This repository contains all the data, python scripts and Jupyter notebooks required to reproduce the analysis and figures presented in the manuscript [The topological shape of gene expression across the evolution of flowering plants](https://www.biorxiv.org/content/10.1101/2022.09.07.506951v1).

We would suggest going through the Jupyter notebooks in the following order:

- `dimension_reduction.ipynb`: Initial data exploration. Creates PCA and t-SNE projections, visualizes them in 2D scatter plots. The scatter plots can be colored by the *stress*, *tissue*, or the *family* labels associated with the samples. NOTE: the last part requires lens function values stored in `saved_lense.csv`.

- `heatmaps.ipynb`: Initial data exploration. Visualizes correlation distance between gene expression profiles as a heatmap. Also creates a hierarchical (average linkage) clustering dendrogram for the samples using the correlation distance as the metric.

- `plant_evo_mapper.ipynb`: Walks you through the process of computing *stress* and *tissue* lenses, constructing mapper graphs, and visualizing them with nodes replaces by pie-charts.
- `plant_evo_mapper_over_cover_grid.ipynb`: Computes mapper graphs over a grid cover parameter values to illustrate the stability of the mapper graph over a range of parameter values.
- `plant_evo_mapper_grouped.ipynb`: Creating and visualizing simplified mapper graphs using manually curated group assignments for nodes of the mapper graph.
- `plant_evo_mapper_lens_correlation_goea.ipynb`: Walks you through the process of performing GO enrichment analysis on subsets of orthogroups that are highly (positively or negatively) correlated with the lens function values.

<hr>

### data
- `/data`: Contains all the inputs, including the gene expression data and metadata files, lens function values, and a color mapping file used to ensure consistent color scheme across all plots, figures and html visualizations.
- `/data/saved_mapper_graphs/*`: contains mapper graphs computed with different lenses and cover parameter values, stored as `.pkl` files.
- `/data/node_group_mappings/*`: contains `*.csv` files which have the mappings between mapper nodes and manually curated group labels, required to create simplified mapper graphs.
- `/data/goea/*`: contains all the files needed to perform GO enrichment analysis.

### code
- `/code/python/`: Contains all the python scripts and Jupyter notebooks.
- `/code/kepler_mapper_pie_nodes/*`: contains `.html`, `.js`, and other necessary files required to visualize mapper graph with nodes replaced by pie charts showing membership distribution for samples in each node.

### results
`/results` and its subdirectories contain all the results: figures, text and `.csv` files produced from the analyses.

<hr>

### Requirements
The code is written for Python 3.8 and has been tested with the following package configuration:

```
goatools==1.2.3
json5==0.9.5
kmapper==2.0.1
matplotlib==3.4.3
numpy==1.20.3
pandas==1.3.3
scikit-learn==1.0
scipy==1.7.1
seaborn==0.11.2
statsmodels==0.12.2
```

Most of these are standard python packages that are well maintained and are backward and forward compatible. So the code should work with older or newer versions of the packages (although we haven't tested it). The only uncommon packages are `kmapper` and `goatools`. For further information and documentation for these packages, check the following links: [kepler mapper](https://kepler-mapper.scikit-tda.org/en/latest/), [goatools](https://pypi.org/project/goatools/).

<hr>

If you find this code useful, please cite our paper:

```
@Article{OstingPalandeWang2020,
	author = {Palande, Sourabh and Kaste, Joshua A.M. and Roberts, Miles D. and Segura Aba, Kenia and Claucherty, Carly and Dacon, Jamell and Doko, Rei and Jayakody, Thilani B. and Jeffery, Hannah R. and Kelly, Nathan and Manousidaki, Andriana and Parks, Hannah M. and Roggenkamp, Emily M. and Schumacher, Ally M. and Yang, Jiaxin and Percival, Sarah and Pardo, Jeremy and Husbands, Aman Y. and Krishnan, Arjun and Montgomery, Beronda L. and Munch, Elizabeth and Thompson, Addie M. and Rougon-Cardoso, Alejandra and Chitwood, Daniel H. and VanBuren, Robert},
	title = {The topological shape of gene expression across the evolution of flowering plants},
	year = {2022},
	doi = {10.1101/2022.09.07.506951},
	publisher = {Cold Spring Harbor Laboratory},
	eprint = {https://www.biorxiv.org/content/early/2022/09/09/2022.09.07.506951.full.pdf},
	journal = {bioRxiv}
```
