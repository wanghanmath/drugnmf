# drugnmf
A tool for drug's adverse effect (AE)-signature extraction and prediction

# In the "code" folder
cmap.R calculates the connectivity score to measure similarity between drugs from expression level.

MDS_imputation.R performs MDS-based network imputation when the drugs' information are imbanlance among different datasets, after network imputation, the dimension of different networks are consistent.

SNF_integration.R applies SNF-based network integration method to integrate multiple dimension-consistent networks.

simulation_blockwise.R and simulation_nonblockwise.R performs simulation under the blockwise and non-blockwise scenarios separately.

# In the "data" folder
"first-stage" folder contains four first-stage networks describing drug-drug similarity from adverse event (AE), bioassay, chemical structure and gene expression level.

"MDS" folder contains four second-stage networks describing drug-drug similarity after MDS-based network imputation, which are the results of MDS_imputation.R 
