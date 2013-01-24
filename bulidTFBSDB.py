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
import time

mgr = Manager()

## Multipocessing function for FIMO
def runFimo(motif, seqFile='footprintSeqs/test.fasta'):
#def runFimo(motif, seqFile='footprintSeqs/footprintSequences.fasta'):
    outFile = open('tmp/pssm_'+motif.replace('$','_')+'.meme','w')
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 3.0\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: + -\n\n'
    memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
    memeHeader += 'A 0.250 C 0.250 G 0.250 T 0.250\n\n'
    outFile.write(memeHeader)
    outFile.write(writeMotifs[motif].getMemeFormatted())
    outFile.close()
    motifHits[motif] = fimo(queryFile='tmp/pssm_'+motif.replace('$','_')+'.meme', seqFile=seqFile)

## Run FIMO on specified files
def fimo(queryFile=None, seqFile=None):
    fimoArgs = '--max-stored-scores 1000000 --verbosity 4 --bgfile bgFile.meme -text --thresh 1e-5 '
    fimoArgs += str(queryFile)+' '+str(seqFile)
    print fimoArgs
    errOut = open('stderr.out','w')
    fimoProc = Popen("~/bin/fimo " + fimoArgs, shell=True,stdout=PIPE,stderr=errOut)
    errOut.close()
    output = fimoProc.communicate()[0].split('\n')
    #print output
    # Write out to a file
    #outFile = open('out.txt','w')
    #outFile.write('\n'.join(output))
    #outFile.close()
    # Return results like cMonkey
    res1 = []
    print output.pop(0)
    #pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence
    for line1 in output:
        splitUp = line1.split('\t')
        #print splitUp
        if len(splitUp)==9:
            res1.append({'chr': splitUp[1], 'start': splitUp[2], 'stop': splitUp[3], 'strand': splitUp[4], 'score':splitUp[5], 'p.value':float(splitUp[6]), 'match.sequence':splitUp[8]})
    return res1

#########################################################
# Load up motifs and turn them into MEME formatted file #
#########################################################
# JASPAR
# pklFile = open('PSSMs/jasparCoreVertebrata_IN_TRANSFAC.pkl','rb')
# jaspar = pklFile.load()
# pklFile.close()

# TransFac
pklFile = open('PSSMs/transfac_2012.1_PSSMs.pkl','rb')
transfac = cPickle.load(pklFile)
pklFile.close()

# UW Novel Motifs
pklFile = open('PSSMs/novelPSSMsUW.pkl','rb')
novel = cPickle.load(pklFile)
pklFile.close()

writeMotifs = mgr.dict(transfac, **novel)
#writeMotifs = dict(novel)
del novel
del transfac

# Prep for run
if not os.path.exists('tmp'):
    os.mkdir('tmp')
motifHits = mgr.dict()

# Run fimo to get the target sites
print 'Running Weeder...'
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.'
pool = Pool(processes=cpus)
#runFimo(writeMotifs.keys()[0])
pool.map(runFimo,writeMotifs.keys()[0:8])


#t1 = time.time()
#res1 = fimo(queryFile='tmp/pssm_'+motif.replace('$','_')+'.meme', seqFile='footprintSeqs/test.fasta')
#res1 = fimo(queryFile='tmp/pssm_'+motif.replace('$','_')+'.meme', seqFile='footprintSeqs/footprintSequences.fasta', numSeqs=float(100))
#t2 = time.time()
#print 'Took %0.3f s' % (t2-t1)
#break

