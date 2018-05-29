###
#
# Federation
# 141105
#
###

# Plot the number of peaks based on the number of samples sequenced

import os
import sys
sys.path.append('/ark/home/af661/src/utils/')
sys.path.append('/ark/home/af661/src/pipeline/')
import pipeline_dfci
import utils
import random

annotationFile = '/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc'
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
projectFolder = '/crusader/projects/cll/final/'

#dataFile = '/ark/home/af661/projects/cll/final/dataTables/CLL_ATAC_TABLE_CLUSTERING.txt'
dataFile = '/crusader/projects/cll/final/dataTables/CLL_H3K27AC_TABLE_NOINPUT.txt'
dataDict = pipeline_dfci.loadDataTable(dataFile)
namesList = dataDict.keys()
namesList = [x for x in dataDict.keys() if x[0] == 'p']
print namesList
numSamples = len(namesList)

matrixDict = {}

#allPeakTable = utils.parseTable('/grail/projects/cll/final/macsEnriched/all_merge.bed', '\t')
#allPeakTable = utils.parseTable('/crusader/projects/cll/final/macsEnriched/pCLL_H3K27AC_merge.bed', '\t')
allPeakTable = utils.parseTable('/crusader/projects/cll/final/venn/merge_beds/CLL_merge.sorted.bed', '\t')
allPeakLoci = [utils.Locus(x[0], x[1], x[2], '.') for x in allPeakTable]
#allPeakCollection = utils.LocusCollection(allPeakLoci, 100)

for name in namesList:

    print name
    matrixDict[name] = {}

    #bedFile = '/crusader/projects/cll/final/macsEnriched/' + name + '_peaks.bed'
    bedFile = '/crusader/projects/cll/final/rose/' + name + '_ROSE/' + name + '_peaks_AllEnhancers.table.txt'
    bed = utils.parseTable(bedFile, '\t')
    #sampleLoci = [utils.Locus(x[0], x[1], x[2], '.') for x in bed]
    sampleLoci = [utils.Locus(x[1], x[2], x[3], '.') for x in bed[6:1006]] #allEnhancer file
    sampleCollection = utils.LocusCollection(sampleLoci, 50)

    for locus in allPeakLoci:
        
        overlap = sampleCollection.getOverlap(locus)
        if overlap:
            matrixDict[name][locus] = 1
        else:
            matrixDict[name][locus] = 0


gainedVEL = []
gainedTable = [['Chr', 'Start', 'End', 'Num_CLL_Samples', 'Num_CD19_Samples']]
for locus in allPeakLoci:
    
    cancerScore = 0
    normalScore = 0

    for name in namesList:
        
        if name[0] == 'p' and matrixDict[name][locus] == 1:
            cancerScore += 1
        if name[0] == 'C' and matrixDict[name][locus] == 1:
            normalScore += 1

#    if cancerScore < 5 and normalScore > 3:
    if True:       
        gainedVEL.append(locus)
        newline = [locus.chr(), locus.start(), locus.end(), cancerScore, normalScore]
        gainedTable.append(newline)

utils.unParseTable(gainedTable, '/crusader/projects/cll/final/vels/sample_count_superEnhancer_top1000.txt', '\t')
