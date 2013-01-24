#################################################################
# @Program: buildTFBSDB.py                                      #
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

###############
### IMPORTS ###
###############
from pssm import pssm
import cPickle, gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
from multiprocessing import Pool, cpu_count, Manager

## Run FIMO on specified files
def fimo(queryFile=None, seqFile=None, numSeqs=None):
    fimoArgs = '--text --verbosity 1 --bgfile bgFile.meme '
    if not numSeqs==None:
        fimoArgs += '--output-pthresh 0.05 ' #+str(float(0.05)/float(numSeqs))+' '
    fimoArgs += str(queryFile)+' '+str(seqFile)
    print fimoArgs
    #errOut = open('tmp/meme/stderr.out','w')
    fimoProc = Popen("~/bin/fimo " + fimoArgs, shell=True,stdout=PIPE) #,stderr=errOut)
    output = fimoProc.communicate()[0].split('\n')
    #print output
    # Write out to a file
    outFile = open('out.txt','w')
    outFile.write('\n'.join(output))
    outFile.close()
    # Return results like cMonkey
    res1 = {}
    #print output.pop(0)
    # Motif	Seq	Start	Stop	Log-odds	p-value	Site
    for line1 in output:
        splitUp = line1.split('\t')
        #print splitUp
        if len(splitUp)==8:
            if not splitUp[1] in res1:
                res1[splitUp[1]] = { str(splitUp[2])+'_'+str(splitUp[3]):{ 'orientation':splitUp[0][0], 'start':splitUp[2], 'stop':splitUp[3], 'strand':splitUp[4], 'logOdds':float(splitUp[5]), 'pValue': float(splitUp[6]), 'site': splitUp[7] } }
            else:
                res1[splitUp[1]][str(splitUp[2])+'_'+str(splitUp[3])] = { 'orientation':splitUp[0][0], 'start':splitUp[2], 'stop':splitUp[3], 'strand':splitUp[4], 'logOdds':float(splitUp[5]), 'pValue': float(splitUp[6]), 'site': splitUp[7] }
    res2 = {}
    for gene in res1:
        for hit in res1[gene]:
            if not gene in res2:
                res2[gene] = res1[gene][hit]['pValue']
            elif res2[gene]>=res1[gene][hit]['pValue']:
                res2[gene] = res1[gene][hit]['pValue']
    #print res2
    if len(res2)>0:
        #res3 = benjaminiHochberg(res2, tests=numSeqs, alpha=0.001)
        res3 = bonferroni(res2, tests=numSeqs, alpha=0.05)
        print len(res2), len(res3)
        return [[i,str(res2[i])] for i in res3]
    else:
        return []


#########################################################
# Load up motifs and turn them into MEME formatted file #
#########################################################
# JASPAR
# pklFile = open('PSSMs/jasparCoreVertebrata_IN_TRANSFAC.pkl','rb')
# jaspar = pklFile.load()
# pklFile.close()

# TransFac
pklFile = open('PSSMs/transfac_2012.1_PSSMs.pkl','rb')
transfac = pklFile.load()
pklFile.close()

# UW Novel Motifs
pklFile = open('PSSMs/novelPSSMsUW.pkl','rb')
novel = pklFile.load()
pklFile.close()

writeMotifs = dict(transfac, **novel)

# Run analysis
if not os.path.exists('tmp'):
    os.mkdir('tmp')
for motif in writeMotifs:
    outFile = open('tmp/pssm_'+motif[0]+'.meme','w')
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 3.0\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: + -\n\n'
    memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
    memeHeader += 'A 0.250 C 0.250 G 0.250 T 0.250\n\n'
    outFile.write(memeHeader)
    outFile.write(transfacPssms[motif[0]].getMemeFormatted())
    outFile.close()
    # Run fimo to get the target sites
    res1 = fimo(queryFile='tmp/pssm_'+motif[0]+'.meme', seqFile='tmp/meme/fasta/eQTL_'+str(motif[1]).replace('_','.')+'.fasta',numSeqs=float(100))
    break

