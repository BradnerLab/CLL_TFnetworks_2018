import os
import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils
import CRC_CLL as crc

from string import upper
from subprocess import call
import numpy
import subprocess

projectFolder = '/ark/home/af661/projects/cll/final/network/'
projectName = 'CLL_NFATC1'
annotationFile = '/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc'
motifConvertFile = '/ark/home/af661/src/coreTFnetwork/annotations/MotifDictionary.txt'
genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
motifDatabaseFile = '/ark/home/af661/src/coreTFnetwork/annotations/VertebratePWMs.txt'

TFfile = '/ark/home/af661/src/coreTFnetwork/annotations/TFlist_NMid_hg19.txt'
TFtable = utils.parseTable(TFfile, '\t')
TFlist = [line[0] for line in TFtable]
TFlistGene = [line[1] for line in TFtable]

refseqToNameDict = {}
annotTable = utils.parseTable(annotationFile, '\t')
for line in annotTable:
    gid = line[1]
    genename = upper(line[12])
    refseqToNameDict[gid] = genename

geneToRefseqDict = {}
for line in annotTable:
    gid = line[1]
    genename = upper(line[12])
    geneToRefseqDict[genename] = gid


# Get the file that contains all enhancer assignments

enhancerFile = '/grail/projects/cll/final/macsEnriched/pCLL_ATAC_merge.bed'
enhancerTable = utils.parseTable(enhancerFile, '\t')
enhancerLoci = [utils.Locus(x[0], x[1], x[2], '.') for x in enhancerTable]

# Gene expression

def createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict, expCutoff, geneToRefseqDict):
    '''
    input: an activity table with refseq in first column and expression or promoter
    acetylation in second column
    output: a dictionary keyed by refseq that points to activity
    '''

    print 'CREATING EXPRESSION DICTIONARY'

    expresionFilename = '/grail/projects/cll/final/expression/CLL_WGS_RNA-Seq.cufflinks.gct'
    expressionTable = utils.parseTable(expresionFilename, '\t')

    expressionDictNM = {}
    expressionDictGene = {}

    for line in expressionTable[3:]:
        geneName = line[1]
        try:
            trid = geneToRefseqDict[geneName]
        except KeyError:
            continue
        values = []
        for sample in line[2:-6]:
            values.append(float(sample))

        exp = numpy.mean(values)

        # Save the expression value of each NMid in a dict, keep higher value if multiple
        if trid in expressionDictNM and exp > expressionDictNM[trid]:
            expressionDictNM[trid] = exp
        elif trid not in expressionDictNM:
            expressionDictNM[trid] = exp

        # Save the value of the expression if it's the highest for that gene
        if geneName in expressionDictGene and exp > expressionDictGene[geneName]:
            expressionDictGene[geneName] = exp
        elif geneName not in expressionDictGene:
            expressionDictGene[geneName] = exp

    cutoff = numpy.percentile(expressionDictGene.values(), expCutoff)
    print 'Expression cutoff: ' + str(cutoff)

    expressedGenes = []
    expressedNM = []

    for nmid in expressionDictNM:
        if float(expressionDictNM[nmid]) > cutoff:
            expressedGenes.append(refseqToNameDict[nmid])
            expressedNM.append(nmid)

    expressedGenes = utils.uniquify(expressedGenes)
    Genefilename = projectFolder + projectName + '_EXPRESSED_GENES.txt'
    utils.unParseTable(expressedGenes, Genefilename, '')

    expressedNM = utils.uniquify(expressedNM)
    NMfilename = projectFolder + projectName + '_EXPRESSED_NM.txt'
    utils.unParseTable(expressedNM, NMfilename, '')

    return expressedNM, expressionDictNM


# Gene assignment


def assignEnhancerToGene(enhancer, tssCollection, startDict, expressionDictNM):

    # If the enhancer overlaps a TSS, save it
    overlappingLoci = tssCollection.getOverlap(enhancer, 'both')
    overlappingGenes =[]
    for overlapLocus in overlappingLoci:
        overlappingGenes.append(overlapLocus.ID())

    # Find all gene TSS within 100 kb
    proximalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,100000,100000),'both')
    proximalGenes =[]
    for proxLocus in proximalLoci:
        proximalGenes.append(proxLocus.ID())

    # If no genes are within 100 kb, find the closest active gene
    closestGene = ''
    if len(overlappingGenes) == 0 and len(proximalGenes) == 0:
        
        distalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,1000000,1000000),'both')
        distalGenes =[]
        for distalLocus in distalLoci:
            distalGenes.append(distalLocus.ID())

        enhancerCenter = (int(enhancer.start()) + int(enhancer.end())) / 2
        distList = [abs(enhancerCenter - startDict[geneID]['start'][0])
                    for geneID in distalGenes]
        if distList:
            closestGene = distalGenes[distList.index(min(distList))]

    overlappingGenes = utils.uniquify(overlappingGenes)
    proximalGenes = utils.uniquify(proximalGenes)
    for refID in overlappingGenes:
        if proximalGenes.count(refID) == 1:
            proximalGenes.remove(refID)


    assignment = ''
    # If a TSS overlaps an enhancer, assign them together
    if overlappingGenes:
        for gene in overlappingGenes:
            if gene in TFlist:
                assignment = gene

    # Otherwise, assign the enhancer to the most active gene in 100 kb
    elif not overlappingGenes and proximalGenes:
        highestGene = ''
        highestActivity = 0
        for gene in proximalGenes:
            if expressionDictNM[gene] > highestActivity:
                highestActivity = expressionDictNM[gene]
                highestGene = gene
        assignment = highestGene

    elif not overlappingGenes and not proximalGenes and closestGene:
        assignment = closestGene

    return assignment

#####
#
# Main
#
#####

# Parse Gene Expression
expCutoff = 50
#expressedNM, expressionDictNM = createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict, expCutoff, geneToRefseqDict)

#startDict = utils.makeStartDict(annotationFile)

#tssLoci = []
#for gene in expressedNM:
#    tssLoci.append(utils.makeTSSLocus(gene,startDict,2500,2500))
#tssCollection = utils.LocusCollection(tssLoci,50)

# Make a FASTA of all ATAC peaks with gene assignments

def makeFASTAfromLoci(enhancerLoci, tssCollection, startDict, genomeDirectory, expressionDictNM):

    fasta = []
    for locus in enhancerLoci:
        
        gene = assignEnhancerToGene(locus, tssCollection, startDict, expressionDictNM)
        
        header = [locus.chr() + '|' + str(locus.start()) + '|' + str(locus.end()) + '|' + gene]
        seq = [utils.fetchSeq(genomeDirectory, locus.chr(), locus.start(), locus.end())]
        
        fasta.append(header)
        fasta.append(seq)

    utils.unParseTable(fasta, '/ark/home/af661/projects/cll/final/fasta/pCLL_ATAC_merge.fasta', '\t')

#makeFASTAfromLoci(enhancerLoci, tssCollection, startDict, genomeDirectory, expressionDictNM)

# Run FIMO

def findMotifs(canidateGenes, projectFolder, projectName, fastaFile, motifConvertFile, motifDatabaseFile):
    '''
    takes the refseq to subpeak seq dict
    returns the networkx object with all connections
    '''

    # Create a dictionary to call motif names keyed on gene names

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]
    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    print 'GENERATING TF NETWORK'

    # select the TF candidates that have motifs
    canidateMotifs = []
    for gene in canidateGenes:
        if gene in motifNames:
            canidateMotifs.append(gene)

    print 'Number of annotated canidate TFs that have motifs: ' + str(len(canidateMotifs))
    canidateMotifs = sorted(canidateMotifs)

    #canidateMotifs = ['NANOG', 'POU5F1', 'SOX2']

    bgCmd = 'fasta-get-markov -m 1 < ' + fastaFile + ' > ' + projectFolder + projectName + '_bg.meme'
    subprocess.call(bgCmd, shell=True)

    utils.formatFolder(projectFolder + 'FIMO/', True)

    fimoCmd = 'fimo'
    for TF in canidateMotifs:
        print TF
        for x in motifDatabaseDict[TF]:
            fimoCmd += ' --motif ' + "'%s'" % (str(x))

    #fimoCmd += ' --thresh 1e-2'
    fimoCmd += ' -verbosity 1'  # thanks for that ;)!
    fimoCmd += ' -text'
    fimoCmd += ' -oc ' + projectFolder + 'FIMO'
    fimoCmd += ' --bgfile ' + projectFolder + projectName + '_bg.meme'
    fimoCmd += ' ' + motifDatabaseFile + ' '
    fimoCmd += fastaFile
    fimoCmd += ' > '+ projectFolder + 'FIMO/fimo.txt'  ##
    print fimoCmd

    fimoOutput = subprocess.call(fimoCmd, shell=True)  #will wait that fimo is done to go on

    return fimoCmd

canidateGenes = ['NFATC1']
fastaFile = '/ark/home/af661/projects/cll/final/fasta/pCLL_ATAC_merge.fasta'
#findMotifs(canidateGenes, projectFolder, projectName, fastaFile, motifConvertFile, motifDatabaseFile)

# Process FIMO

fimo = utils.parseTable('/grail/projects/cll/final/network/FIMO/fimo.txt', '\t')

motifDatabase = utils.parseTable(motifConvertFile, '\t')
motifDatabaseDict = {}
motifNames = [line[1] for line in motifDatabase]
for line in motifDatabase:
    motifDatabaseDict[line[0]] = line[1]

motifCount = {}
for line in fimo[1:]:

    source = motifDatabaseDict[line[0]]   #motifId
    region = line[1].split('|')

    if region[3]:
        target = refseqToNameDict[region[3]]   #gene name corresponding to the NMid
    else:
        target = False


    motifLocus = utils.Locus(region[0], int(region[1]) + int(line[2]), int(region[1]) + int(line[3]), '.')

    if target and target not in motifCount:
        motifCount[target] = [motifLocus]
        
    elif target:

        # check if motif is in a unique location
        collection = utils.LocusCollection(motifCount[target])
        if not collection.getOverlap(motifLocus):
              motifCount[target].append(motifLocus)

    
countTable = []
for gene in motifCount:

    newline = [gene, len(motifCount[gene])]
    countTable.append(newline)

utils.unParseTable(countTable, '/grail/projects/cll/final/network/CLL_NFATC1_motifCount.txt', '\t')
