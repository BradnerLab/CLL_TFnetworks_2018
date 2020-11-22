#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: extendedSuperNetwork.py 
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)
***********************************************
"""
def main():
    # First, load in the node TFs, ATAC peaks and super enhancer regions we'll consider for this analysis
    # From networks already constructed from CRC2.py

    # node_file = '/crusader/projects/cll/final/network/lines/zinba/' + projectName + '/' + projectName + '_NODELIST.txt'
    node_table = utils.parseTable(node_file, '\t')
    nodelist = [x[0] for x in node_table]
    print(nodelist)
    # super_enhancer_file = '/crusader/projects/cll/final/rose/' + projectName + '_H3K27ac/' + projectName + '_H3K27ac_peaks_SuperEnhancers.table.txt'

    se_table = utils.parseTable(super_enhancer_file, '\t')

    # subpeak_file = '/crusader/projects/cll/final/zinba/lines/MEC1_ATAC/MEC1_ATAC.peaks.bed'
    subpeak_table = utils.parseTable(subpeak_file, '\t')
    subpeak_loci = []
    for line in subpeak_table:
        subpeak_loci.append(utils.Locus(line[0], line[1], line[2], '.'))
    subpeak_collection = utils.LocusCollection(subpeak_loci, 100)
    subpeak_dict = {} # key is enhancer ID, points to a list of loci

    # assign subpeak Loci to each super enhancer
    fasta = []
    se_namelist = []
    for line in se_table[6:]:

        se_id = line[0]
        se_namelist.append(se_id)
        subpeak_dict[se_id] = []
        
        se_locus = utils.Locus(line[1], line[2], line[3], '.')
        overlaps = subpeak_collection.getOverlap(se_locus)

        for overlap in overlaps:
            subpeak_dict[se_id].append(overlap)

            subpeak = overlap
            
            fastaTitle = se_id + '|'  + subpeak.chr() + '|' + str(subpeak.start()) + '|' + str(subpeak.end())
            fastaLine = utils.fetchSeq(genomeDirectory, subpeak.chr(), int(subpeak.start()+1), int(subpeak.end()+1))

            fasta.append('>' + fastaTitle)
            fasta.append(upper(fastaLine))


    subpeaksFasta = "{0}/{1}_SUBPEAKS.fa".format(projectFolder,projectName)
    utils.unParseTable(fasta, subpeaksFasta, '')


    # call FIMO and find the motifs within each enhancer

    # motifConvertFile = '/ark/home/af661/src/coreTFnetwork/annotations/MotifDictionary.txt'
    # motifDatabaseFile = '/ark/home/af661/src/coreTFnetwork/annotations/VertebratePWMs.txt'

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]
    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    canidateMotifs = []
    for gene in nodelist:
        if gene in motifNames:
            canidateMotifs.append(gene)
    print(canidateMotifs)

    # bgCmd = 'fasta-get-markov -m 1 < ' + projectFolder + projectName + '_SUBPEAKS.fa > ' + projectFolder + projectName + '_bg.meme'
    bgMeme = "{0}/{1}_bg.meme".format(projectFolder, projectName)
    bgCmd = "fasta-get-markov -m 1 < {0} > {1}".format(subpeaksFasta, bgMeme)
    subprocess.call(bgCmd, shell=True)

    fimoFolder = "{0}/FIMO".format(projectFolder)
    fimoFile  =  "{0}/fimo.txt".format(fimoFolder)
    utils.formatFolder(fimoFolder, True)

    fimoCmd = 'fimo'
    for TF in canidateMotifs:
        print(TF)
        for x in motifDatabaseDict[TF]:
            # fimoCmd += ' --motif ' + "'%s'" % (str(x))
            fimoCmd = "{0} --motif {1}".format(fimoCmd, x)

    # ##fimoCmd += ' --thresh 1e-5'
    # fimoCmd += ' -verbosity 1'  # thanks for that ;)!
    # fimoCmd += ' -text'
    # fimoCmd += ' -oc ' + projectFolder + 'FIMO'
    # fimoCmd += ' --bgfile ' + projectFolder + projectName + '_bg.meme'
    # fimoCmd += ' ' + motifDatabaseFile + ' '
    # fimoCmd += projectFolder + projectName + '_SUBPEAKS.fa'
    # fimoCmd += ' > '+ projectFolder + 'FIMO/fimo.txt'  ##
    fimoCmd = "{0} -verbosity 1 -text -oc {1} --bgfile {2} {3} {4} > {5}".format(fimoCmd, fimoFolder, bgMeme, motifDatabaseFile, subpeaksFasta, fimoFile)
    print(fimoCmd)

    fimoOutput = subprocess.call(fimoCmd, shell=True)  #will wait that fimo is done to go on


    # next, build a dictionary with all network info and output a matrix with the same information
    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]

    # The reverse of the other dict, from motif name to gene name
    for line in motifDatabase:
        motifDatabaseDict[line[0]] = line[1]


    fimoTable = utils.parseTable(fimoFile, '\t')
    motifDict = defaultdict(dict)
    for line in fimoTable[1:]:
        source = motifDatabaseDict[line[0]]   #motifId
        region = line[1].split('|')
        target = region[0]   #gene name corresponding to the NMid

        if target not in motifDict[source]:
            motifDict[source][target] = 0
        motifDict[source][target] += 1

    # make matrix
    matrix = [se_namelist]
    for tf in motifDict:
        newline = [tf]
        for se_id in se_namelist:
            if se_id in motifDict[tf]:
                newline.append(motifDict[tf][se_id])
            else:
                newline.append(0)
        matrix.append(newline)

    matrix_name = "{0}/{1}_extendedNetwork.allEnhancers.matrix.txt".format(projectFolder, projectName)
    utils.unParseTable(matrix, matrix_name, '\t')

################ USER DEFINED FUNCTIONS ###################
def print_help():
    ''' Print system help '''
    print >> sys.stderr, "\n ----------------- HELP ------------------\n", parser.print_help(), "\n"

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python extendedSuperNetwork.py -pn=127abcam-5320-PPT-1 -ad=/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/analysis/crcs/extendedNetworks -nf=/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/analysis/crcs/127abcam-5320-PPT-1/crcs127abcam-5320-PPT-1_NODELIST.txt -se=/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/analysis/crcs/127abcam-5320-PPT-1/rose2/127abcam-5320-PPT-1_summits_SuperEnhancers.table.txt -ap=/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/analysis/customPeaks/mergedLibrary/q0_01/5320_PPT-1_004_R1_peaks.broadPeak -gn=MM10
        -------------------------------------------------
        CONTACT: 
            Gaurav Jain
            gaurav.jain@tum.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-pn", metavar='--prname', help="*Project name"    , dest="projectName"  , type=str, required=True)
    parser.add_argument("-ad", metavar='--anydir', help="*Analysis directry", dest="analysisDir", type=str, required=True)
    parser.add_argument("-nf", metavar='--nodefl', help="*Node file from the CRC analysis. It's <filename>_NODELIST.txt", dest="node_file", type=str, required=True)
    parser.add_argument("-se", metavar='--supenf', help="*Super enhancers file from the CRC analysis. It's crcdir/rose2/<filename>_SuperEnhancers.table.txt", dest="super_enhancer_file", type=str, required=True)
    parser.add_argument("-ap", metavar='--atacpk', help="*Atacseq peaks file. It is a BED file of regions to search for motifs", dest="subpeak_file", type=str, required=True)
    parser.add_argument("-mt", metavar="--motifs", help=" Enter an alternative PWM file for the analysis", dest="motifs")
    parser.add_argument("-gn", metavar="--genome", help=" Provide the build of the genome to be used for the analysis. Currently supports HG19, HG18, MM9 and MM10\n Default = MM10", dest="genome", default="MM10")
    parser.add_argument("-uc",         "--ucscfm", help=" If set, use the ucsc folders or files with chromosome names as chr1, chr2, etc.", dest="is_ucsc", action='store_true', default=False,)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().node_file:
        logdir="{0}/logs".format(get_file_info(parser.parse_args().analysisDir)[0])
        create_dir(logdir)
        logfile = "{0}/{1}_extended_super_network.log".format(logdir, get_file_info(parser.parse_args().node_file)[1])
        print(logdir)
    else:
        logdir  = "{0}/logs".format(os.getcwd())
        create_dir(logdir)
        logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])

    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args
    
if __name__=="__main__":
    print (__doc__)

    # Built in modules
    import os
    import sys
    import argparse
    sys.path.insert(0,'/home/rad/users/gaurav/projects/ctrc/scripts/pipeline')
    import utils

    import string
    import textwrap
    import numpy as np
    import scipy as sp
    import scipy.stats as stats

    import subprocess
    import os

    # from string import upper
    from random import randrange
    from collections import defaultdict

    import networkx as nx
    from networkx.algorithms.clique import find_cliques_recursive
    import pickle

    # user defined modules
    from gjainPyLib import *      # import all the functions from the Gaurav`s python library


    ################ USER CONFIGURATION ###################
    np.set_printoptions(precision=6)
    #######################################################

    # Get input options
    args = check_options()

    # Store the variables
    projectName         = args.projectName
    analysisDir       = args.analysisDir
    node_file           = args.node_file
    super_enhancer_file = args.super_enhancer_file
    subpeak_file        = args.subpeak_file
    motifs              = args.motifs
    genome              = args.genome
    is_ucsc             = args.is_ucsc
    # genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
    genome = genome.upper()
    if genome == 'HG19':
        genomeDirectory = '/home/rad/packages/data/fasta/human/hg19/chromosomes/'
        annotationFile = '/home/rad/users/gaurav/projects/ctrc/scripts/pipeline/annotation/hg19_refseq.ucsc'
        TFfile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/TFlist_NMid_hg19.txt'

    if genome == 'HG18':
        genomeDirectory = '/grail/genomes/Homo_sapiens/human_gp_mar_06_no_random/fasta/'
        annotationFile = '/ark/home/cl512/src/pipeline/annotation/hg18_refseq.ucsc'
        TFfile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/TFlist_NMid_hg19.txt'

    if genome == 'MM9':
        genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
        annotationFile = '/home/rad/users/gaurav/projects/ctrc/scripts/pipeline/annotation/mm9_refseq.ucsc'
        TFfile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/TFlist_NMid_mm9.txt'

    if genome == 'MM10':
        TFfile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/TFlist_NMid_mm10.txt'
        if is_ucsc:
            genomeDirectory = '/home/rad/packages/data/fasta/mouse/mm10/ucsc_chromosomes/'
            annotationFile = '/home/rad/users/gaurav/projects/ctrc/scripts/pipeline/annotation/ucsc/mm10_refseq.ucsc'
        else:
            genomeDirectory = '/home/rad/packages/data/fasta/mouse/mm10/chromosomes/'
            annotationFile = '/home/rad/users/gaurav/projects/ctrc/scripts/pipeline/annotation/mm10_refseq.ucsc'

    if motifs:
        motifDatabaseFile = options.motifs
    else:
        motifConvertFile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/MotifDictionary.txt'
        motifDatabaseFile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/VertebratePWMs.txt'

    # projectFolder = '/crusader/projects/cll/final/network/extended/node-SE/'
    # projectName = 'MEC1'
    # utils.formatFolder(projectFolder + projectName, True)
    # projectFolder = projectFolder + projectName + '/'
    projectFolder = "{0}/{1}".format(analysisDir, projectName); create_dir(projectFolder)
    main()


