import os
import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

import string

import numpy
import scipy
import scipy.stats

import subprocess
import os

from string import upper
from random import randrange
from collections import defaultdict

import networkx as nx
from networkx.algorithms.clique import find_cliques_recursive
import pickle

genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'

projectFolder = '/crusader/projects/cll/final/network/extended/node-SE/'
projectName = 'MEC1'
utils.formatFolder(projectFolder + projectName, True)
projectFolder = projectFolder + projectName + '/'

# First, load in the node TFs, ATAC peaks and super enhancer regions we'll consider for this analysis
# From networks already constructed from CRC2.py
nodelist = ['PAX5', 'SMAD3', 'IKZF1', 'RUNX3', 'TBX21', 'NR3C1', 'FOXO1', 'FLI1', 'RARA', 'RREB1', 'KLF13', 'STAT1', 'BHLHE40', 'TCF3', 'IRF2', 'ETV6', 'MAX', 'HIVEP2', 'ELF1', 'IKZF2', 'IRF4']

#super_enhancer_file = '/crusader/projects/cll/final/rose/pCLL_' + projectName + '_H3K27ac_ROSE/pCLL_' + projectName + '_H3K27ac_peaks_SuperEnhancers.table.txt'
#super_enhancer_file = '/crusader/projects/cll/final/rose/pCLL_' + projectName + '_H3K27ac_ROSE/pCLL_' + projectName + '_H3K27ac_peaks_SuperEnhancers.table.txt'
super_enhancer_file = '/crusader/projects/cll/final/rose/' + projectName + '_H3K27ac/' + projectName + '_H3K27ac_peaks_SuperEnhancers.table.txt'
se_table = utils.parseTable(super_enhancer_file, '\t')

matrix_file = projectFolder + projectName + '_extendedNetwork.matrix.txt'
matrix_table = utils.parseTable(matrix_file, '\t')

gene_target_file = '/crusader/projects/cll/final/rose/' + projectName + '_H3K27ac/' + projectName + '_H3K27ac_peaks_SuperEnhancers_ENHANCER_TO_TOP_GENE.txt'
gene_target_table = utils.parseTable(gene_target_file, '\t')
gene_target_dict = {}
for x in gene_target_table[1:]:
    gene_target_dict[x[0]] = x[12]
print  gene_target_dict

edict = {}
enh_list = matrix_table[0]
for enh in matrix_table[0]:
    edict[enh] = {}

for line in matrix_table[1:]:
    tf = line[0]
    for i, enh in enumerate(enh_list):
        edict[enh][tf] = line[i+1]


#print edict['6_MACS_peak_8891_lociStitched']['ZNF217']

output = [['# Listing enhancers bound by these genes:'], nodelist, ['ID', 'TopGene']]
for enh in enh_list:
    score = 0
    for tf in nodelist:
        if int(edict[enh][tf]) > 0:
            score += 1

    if score == len(nodelist):

        newline = [enh, gene_target_dict[enh]]
        output.append(newline)

utils.unParseTable(output, projectFolder + projectName + '_CRC-VSA_targetsBCL2.txt', '\t')
