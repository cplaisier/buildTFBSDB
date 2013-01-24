#################################################################
# @Program: compileSequences.py                                 #
# @Version: 1                                                   #
# @Author: Christopher L Plaisier, PhD                          #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
# American Cancer Society Postdoctoral Fellowship               #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  1/15/2013                      #
#################################################################

import os, gzip
from copy import deepcopy
import tarfile

################################################
# 1. Read in footprints (8,374,968 footprints) #
################################################
entries = 0
footprints = {} # chr -> [footprint, ... ] - ordered by genomic region
inFile = open('combined.fps','r')
# Read in line by line and build footprints dictionary
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.strip().split()
    if not splitUp[0] in footprints:
        footprints[splitUp[0]] = []
    footprints[splitUp[0]].append([splitUp[0], int(splitUp[1]), int(splitUp[2])])
    entries += 1
inFile.close()
print 'Footprints =',entries

## Merge overlapping sequences 6bp
# For each chromosome
entMerg = 0
fpMerged ={}
for i in footprints:
    fpMerged[i] = []
    prev = ''
    fpLen = len(footprints[i])
    cur = 0
    while 1:
        cur += 1
        if not cur<fpLen:
            break
        if not prev=='' and prev[2]<=(footprints[i][cur][1]-6):
            prev[2] = footprints[i][cur][2]
        else:
            if not prev=='':
                fpMerged[i].append(prev)
                entMerg += 1
            prev = footprints[i][cur]

print 'Footprints =',entries,'; Merged = ',entMerg,'; Merge Length = 6'

##########################
# 2. Parse out sequences #
##########################
print '  Extracting the sequence data...'
# Unzip sequences for extraction
#tar = tarfile.open('sequences/fasta/chromFa.tar.gz')
#tar.extractall(path='sequences/fasta')
#tar.close()

# 6. Extract the sequences
for chrom in fpMerged:
    footprintSeqFile = open('footprintSequences_'+chrom+'.fasta','w')
    if os.path.exists('sequences/fasta/'+str(chrom)+'.fa'):
        chrSeqFile = open('sequences/fasta/'+str(chrom)+'.fa','r')
    elif os.path.exists('sequences/fasta/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa'):
            chrSeqFile = open('sequences/fasta/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa','r')
    else:
        print 'FATAL ERROR!!!! Arghhh',chrom,'(',str(chrom).lstrip('chr').replace('_random',''),')does not have a seqeunce file!'
        break
    chrSeqFile.readline() # Get rid of header
    chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
    chrSeq = ''.join(chrSeq)
    print '  ',chrom,' (',len(chrSeq),')...'
    for fp1 in fpMerged[chrom]:
        footprintSeq = chrSeq[(fp1[1]-1):(fp1[2]-1)]
        footprintSeqFile.write('>'+str(fp1[0])+':'+str(fp1[1])+'-'+str(fp1[2])+'\n'+str(footprintSeq)+'\n')
    footprintSeqFile.close()

