#####
#
# Looking at the effect of JQ1 on nascent RNA-seq at super enhancers
#
#####


import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

from collections import defaultdict
from string import upper
import numpy as np
from math import log

# Annotation file for hg19

annotationFile = '/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc'
startDict = utils.makeStartDict(annotationFile)

print 'making TSS loci'
tssLoci = []
counter = 0
for gene in startDict:
    counter += 1
    if counter%1000 == 0:
        print counter
    tssLoci.append(utils.makeTSSLocus(gene,startDict,100000,100000)) # proximal = within 100kb
tssCollection = utils.LocusCollection(tssLoci,200)

print 'converting gene names'
refseqToNameDict= {}
annotTable = utils.parseTable(annotationFile, '\t')
for line in annotTable:
    gid = line[1]
    genename = upper(line[12])
    refseqToNameDict[gid] = genename

nameToRefseqDict = defaultdict(list)
annotTable = utils.parseTable(annotationFile, '\t')
for line in annotTable:
    gid = line[1]
    genename = upper(line[12])
    nameToRefseqDict[genename].append(gid)

# Define cell type and enhancer files

cell_type = 'OSU'

superFile = '/crusader/projects/cll/final/rose/OSUCLL_H3K27ac/OSUCLL_H3K27ac_peaks_SuperEnhancers_ENHANCER_TO_GENE.txt'
superTable = utils.parseTable(superFile, '\t')
superGenes = []
for line in superTable:
    names = upper(line[11]).split(',') # 11 is proximal genes
    names += upper(line[12])
    for n in names:
        if n not in superGenes:
            superGenes.append(n)

print superGenes
                         
typicalGenes = []

allEnhancersFile = '/crusader/projects/cll/final/rose/OSUCLL_H3K27ac/OSUCLL_H3K27ac_peaks_AllEnhancers.table.txt'
allTable = utils.parseTable(allEnhancersFile, '\t')
allLoci = [utils.Locus(x[1], x[2], x[3], '.') for x in allTable[7:]]

for locus in allLoci:
    overlap = tssCollection.getOverlap(locus)
    if overlap:
        for o in overlap:
            typical_name = refseqToNameDict[o.ID()]
            if typical_name not in superGenes and typical_name not in typicalGenes:
                typicalGenes.append(typical_name)

print typicalGenes
print len(typicalGenes)
allGenes = superGenes + typicalGenes

# Pull in expression data

all_fold =  []
typical_fold = []
super_fold = []

expression_cutoff = 5

expression_table = utils.parseTable('/crusader/projects/cll/final/se_vs_te/CLL_JQ1_treat_all_fpkm_exprs_norm.txt', '\t')

expression_dict = {}


TF_table = utils.parseTable('/crusader/projects/cll/final/se_vs_te/MEC1_RSA_parsed.txt', '\t')
TF_list = [x[0] for x in TF_table]

outtable = []

for line in expression_table[1:]:
    
    gene = line[0]

    newline = [gene]
        
    dmso_values = [float(i) for i in line[13:16]]
    jq1_values = [float(i) for i in line[16:19]]
    
    dmso_mean = np.mean(dmso_values)
    jq1_mean = np.mean(jq1_values)

    newline.append(dmso_mean)
    newline.append(jq1_mean)
    outtable.append(newline)
    if dmso_mean > expression_cutoff or jq1_mean > expression_cutoff:

        fold = float(jq1_mean) / dmso_mean
        log2fold = log(fold, 2)
        
        all_fold.append(fold)
            
        if gene in typicalGenes:
            typical_fold.append(fold)
        if gene in superGenes:
            super_fold.append(fold)


            

outfilename = '/crusader/projects/cll/final/expression/jq1/' + cell_type + 'enhancer_type_foldChange_nodelist.txt'

outfile = []
outfile.append(all_fold)
outfile.append(typical_fold)
outfile.append(super_fold)

utils.unParseTable(outfile, outfilename, '\t')
