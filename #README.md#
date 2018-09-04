# CLL_TFnetworks_2018

## Overview

This document contains scripts used to analyze data and generate figures for Ott et al (2018) in Cancer Cell. For inquries about specific analyses, please contact Alex Federation at federation@altius.org or Christopher Ott at christopher.ott@mgh.harvard.edu. In this README file, we provide either a detailed description  or the exact scripts used to generate each figure. For an up-to-date distribution of the code for network construction, please use pip (https://pypi.org/project/coltron/). For additional information about supplemental figures, please contact Alex Federation.

## Figure 1

### A

Super enhancer plots are generated with the Rank-Ordering of Super Enhancers algorithm, version 2.0. An up-to-date version is found here: https://github.com/BradnerLab/pipeline/ROSE2_main.py

### B

Generated with teh UCSC genome browser

### C

Peaks were called with MACS, as outlined in the methods. Then peaks from each sample were merged into a master dataset and overlaps were calculated using BEDOPS

## Figure 2

### A

Values for input into the heatmap were generated with diffSupers.R using super enhancer tables output from ROSE2.0.

### B, C

Gene tracks were generated with BamplotTurbo. An up-to-date version is found here: https://github.com/BradnerLab/pipeline/bamPlot_turbo.py

### D

Plots generated with diffSeRna.R using the output from diffSupers.R

## Figure 3

### B

Plots are generated as a standard output of the network construction algorithm contained in COLTRON.

### C

Target genes of CRCs were assigned using extendedSuperNetwork.py. Input files include super enhancer data from ROSE2.0, ATAC peak files called from Zinba and annotation files downloaded from the UCSC browser or motif databases. The output is a large matrix with the number of motif occurances for each node TF in the CLL sample within ATAC peaks found in each super enhancers present in that sample. Enhenacers were assiged to genes using targetGenes.py

### D

Clique scores are generated as a standard output of the	network	construction algorithm contained in COLTRON, with the file name *_ENRICHED_CLIQUE_FACTORS.txt

## Figure 4

### A

Clique scores are generated as a standard output of the network construction algorithm contained in COLTRON, with the file name	*_ENRICHED_CLIQUE_FACTORS.txt

### B

RSA p values were calculated using previously published code (König et al., 2007)

## Figure 5

### A

Boxplots were generated using JQ1_RNAseq.py, using annotation files, ROSE2.0 output and normalized RNA-seq FPKM expression data.