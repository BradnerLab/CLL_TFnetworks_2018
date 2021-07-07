import sys
sys.path.insert(0,'/home/rad/users/gaurav/projects/ctrc/scripts/pipeline')
# sys.path.append('/ark/home/af661/src/utils/')
import utils

from string import upper

TFfile = '/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations/TFlist_NMid_hg19.txt'
TFtable = utils.parseTable(TFfile, '\t')
TFlist = utils.uniquify([line[1] for line in TFtable])

annotationFile = '/home/rad/users/gaurav/projects/ctrc/scripts/pipeline/annotation/mm10_refseq.ucsc'
annotTable = utils.parseTable(annotationFile, '\t')

output = []

for line in annotTable:
    gid = line[1]
    genename = upper(line[12])

    if genename in TFlist:
        output.append([gid, genename])

utils.unParseTable('/home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018/annotations','TFlist_NMid_mm10.txt', '\t')