import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

from string import upper

TFfile = '/ark/home/af661/src/coreTFnetwork/annotations/TFlist_NMid_hg19.txt'
TFtable = utils.parseTable(TFfile, '\t')
TFlist = utils.uniquify([line[1] for line in TFtable])

annotationFile = '/ark/home/cl512/pipeline/annotation/mm9_refseq.ucsc'
annotTable = utils.parseTable(annotationFile, '\t')

output = []

for line in annotTable:
    gid = line[1]
    genename = upper(line[12])

    if genename in TFlist:
        output.append([gid, genename])

utils.unParseTable(output, 'TFlist_NMid_mm9.txt', '\t')
