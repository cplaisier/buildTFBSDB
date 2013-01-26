#################################################################
# @Program: convertSELEXMotifs.py                               #
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
# Copyrighted by Chris Plaisier  1/24/2013                      #
#################################################################

###############
### IMPORTS ###
###############
import cPickle, os, re, math
from pssm import pssm
from copy import deepcopy
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *

# Make the files for a TomTom run
def makeQueryFile(nucFreqs, queryPssms, num, strands='+ -'):
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 4\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: '+strands+'\n\n'
    memeHeader += 'Background letter frequencies\n'
    memeHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))
    # Make query PSSM file
    queryFile = open('tmp/query'+str(num)+'.tomtom','w')
    queryFile.write(memeHeader)
    queryFile.write('\n\n'+'\n\n'.join([pssm1.getMeme4Formatted() for pssm1 in queryPssms]))
    queryFile.close()

# Make the files for a TomTom run
def makeTargetFile(nucFreqs, targetPssms, num, strands='+ -'):
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 4\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: '+strands+'\n\n'
    memeHeader += 'Background letter frequencies\n'
    memeHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))
    # Make target PSSM file
    targetFile = open('tmp/target'+str(num)+'.tomtom','w')
    targetFile.write(memeHeader)
    targetFile.write('\n\n'+'\n\n'.join([pssm1.getMeme4Formatted() for pssm1 in targetPssms]))
    targetFile.close()

# Run TomTom on the files
def TomTom(num, distMeth='ed', qThresh='0.05', minOverlap=6):
    # Arguments for tomtom
    tomtomArgs = ' -dist '+str(distMeth)+' -o tmp/tomtom_out -text -thresh '+str(qThresh)+' -min-overlap '+str(minOverlap)+' tmp/query'+str(num)+'.tomtom tmp/target'+str(num)+'.tomtom'
    #print tomtomArgs
    #p = Popen("tomtom" + tomtomArgs, shell=True)
    #sts = os.waitpid(p.pid, 0)
    errOut = open('tmp/stderr.out','w')
    tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE, stderr=errOut)
    outputFile = open('tmp/tomtom/tomtom'+str(num)+'.out','w')
    output = tomtomProc.communicate()[0]
    outputFile.write(output)
    outputFile.close()
    errOut.close()

# Wrapper function to run TomTom using multiprocessing pool
def runTomTom(i):
    TomTom(i, distMeth='ed', qThresh='0.001', minOverlap=6) #blic5

# Method returning the information content of a motif.
def ic(pssm,norm=True):
    bgFreq = [0.21, 0.29, 0.29, 0.21]
    res=0
    pwm=pssm.getMatrix()
    for i in range(len(pwm)):
        res+=2
        for a in [0,1,2,3]:
            if pwm[i][a]!=0:
                if norm==True:
                    res+=pwm[i][a]*math.log(pwm[i][a]/bgFreq[a],2)
                else:
                    res+=pwm[i][a]*math.log(pwm[i][a],2)
    return res

# 1. Read in de novo motif detected PWMs from Uniprobe Protein Binding Microarray (PBM) experiments
uniprobePSSMs = {}
genes = {}
geneNum = []
for dir1 in os.listdir('All_PWMs'):
    for file1 in os.listdir('All_PWMs/'+dir1):
        inFile = open('All_PWMs/'+dir1+'/'+file1,'r')
        inLines = [line for line in inFile.readlines() if line.strip()]
        gene = file1.replace('_pwm','').replace('.pwm','').replace('.txt','').replace('_primary','').replace('_secondary','').replace('.1','').replace('.2','')
        splitUp = gene.split('_')
        if len(splitUp)>1:
            gene = splitUp[0]
        # If enrty already in there increment the .* by one
        if gene in genes:
            last = genes[gene][-1]
            lastNum = str(int(last.split('.')[1])+1)
            genes[gene].append(gene+'.'+lastNum)
            gene = gene+'.'+lastNum
        else:
            genes[gene] = [gene+'.1']
            gene = gene+'.1'
        # Get rid of the header crap
        for i in range(len(inLines)-4):
            crap = inLines.pop(0)
        tmpPssm = []
        for line in inLines:
            splitUp = [i for i in line.strip().split('\t') if i]
            catchMe = splitUp.pop(0)
            for j in range(len(splitUp)):
                if catchMe=='A:':
                    tmpPssm.append([float(splitUp[j])])
                else:
                    tmpPssm[j].append(float(splitUp[j]))
        uniprobePSSMs[gene] = pssm(biclusterName=gene, pssm=deepcopy(tmpPssm), nsites=str(100), eValue=str(0.01))
print 'PSSMs recovered:', len(uniprobePSSMs)
print 'Genes recovered:',len(genes)
red1 = 0
for gene in genes:
    if len(genes[gene])>1:
        red1 += 1
print 'Redundant genes:',red1

# Test to show that the information content produces the expected values
# The code for this funciton was obtained from the Biopyhton Bio.motif.ic funciton
#tmpPssm = [[0.05, 0.05, 0.85, 0.05],
#           [0.85, 0.05, 0.05, 0.05],
#           [0.05, 0.05, 0.85, 0.05],
#           [0.65, 0.05, 0.25, 0.05],
#           [0.85, 0.05, 0.05, 0.05]]
#testPssm = pssm(biclusterName='test', pssm=deepcopy(tmpPssm), nsites=str(100), eValue=str(0.01))
#print 'ic.norm = ',ic(testPssm,norm=True)
#print 'ic = ', ic(testPssm,norm=False) # Biopython Bio.moitf.ic gives 5.27 as the IC

# 3. Write out pickle of the PSSMs for the Uniprobe Protein Binding microarray set
outFile = open('uniprobePSSMs.pkl','wb')
cPickle.dump(uniprobePSSMs,outFile)
outFile.close()

# 4. Compare motifs to identify non-redundant set
if not os.path.exists('tmp/tomtom'):
    os.makedirs('tmp/tomtom')
print 'Making files...'
tested = []
geneNames = genes.keys()
for i in range(len(genes)):
    tmp = genes[geneNames[i]]
    if not os.path.exists('tmp/query'+str(i)+'.tomtom'):
        makeQueryFile(nucFreqs={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, queryPssms=[uniprobePSSMs[tmp.pop(0)]],num=i)
    if not os.path.exists('tmp/target'+str(i)+'.tomtom'):
        makeTargetFile(nucFreqs={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, targetPssms=[uniprobePSSMs[pssm] for pssm in tmp], num=i)
print 'Done.'

# Run this using all cores available
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.'
#runTomTom(0)
pool = Pool(processes=cpus)
pool.map(runTomTom,range(len(genes)))
print 'Done with Tomtom runs.\n'


print 'Removing redundancy in PSSMs...'
redundantPSSMs = []
for run in range(len(genes)):
    # Open file
    outputFile = open('tmp/tomtom/tomtom'+str(run)+'.out','r')
    output = outputFile.readlines()
    outputFile.close()
    # Now iterate through output and save data
    output.pop(0) # Get rid of header
    #Query ID	Target ID	Optimal offset	p-value	E-value	q-value	Overlap	Query consensus	Target consensus	Orientation
    if len(output)>0:
        tmpPssms = []
        # Get ids for redundant motifs
        for i in range(len(output)):
            splitUp = output[i].strip().split('\t')
            # P-value cutoff
            if float(splitUp[3])<=0.00001:
                if len(tmpPssms)==0:
                    tmpPssms.append(splitUp[0])
                tmpPssms.append(splitUp[1])
        # Decision tree for choosing non-redundant examplar
        # 1. Choose the one with the highest information content
        top = ''
        topIc = ''
        for pssm1 in tmpPssms:
            testIc = ic(uniprobePSSMs[pssm1],norm=True)
            if top=='' or testIc>topIc:
                top = pssm1
                topIc = testIc
        redundantPSSMs += [i for i in tmpPssms if not i==top]
        exemplar = top

# Now clean out redundant PSSMs
for pssm1 in redundantPSSMs:
    del uniprobePSSMs[pssm1]
print 'Non-redundant PSSMs: ',len(uniprobePSSMs)
pklFile = open('uniprobePSSMsNonRedundant.pkl','wb')
cPickle.dump(uniprobePSSMs, pklFile)
pklFile.close()

