import cPickle, os, re
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
    print tomtomArgs
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

# 1. Read in de novo motif detected PWMs from systematic SELEX experiments
inFile = open('selex.csv','r')
inFile.readline() # Get rid of header line
# symbol,family,clone type,ligand sequence,batch,Seed,multinomial,cycle,site type,comment,Matrix is  one of the representative PWMs
inFile.readline() # Get rid of header line
selexPSSMs = {}
pssms = {}
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    name = splitUp[0]+'_'+splitUp[1]+'_'+splitUp[2]+'_'+splitUp[8]+'_'+str(len(splitUp[5]))
    if not name in pssms:
        pssms[name] = 1
    else:
        pssms[name] += 1
    name = name+'_'+str(pssms[name])
    # Read in PSSM
    tmpPssm = []
    for cur in ['A','C','G','T']:
        line = inFile.readline()
        splitUp = [i for i in line.strip().split(',') if i]
        catchMe = splitUp.pop(0)
        for j in range(len(splitUp)):
            if cur=='A':
                tmpPssm.append([splitUp[j]])
            else:
                tmpPssm[j].append(splitUp[j])
    # Convert counts to frequencies
    for i in range(len(tmpPssm)):
        tmpSum = float(tmpPssm[i][0])+float(tmpPssm[i][1])+float(tmpPssm[i][2])+float(tmpPssm[i][3])
        for j in [0,1,2,3]:
            tmpPssm[i][j] = float(tmpPssm[i][j])/tmpSum
    # Instantiate PSSM object
    selexPSSMs[name] = pssm(biclusterName=name, pssm=deepcopy(tmpPssm), nsites=str(100), eValue=str(0.01))
inFile.close()
print 'PSSMs recovered:',len(selexPSSMs)

# 3. Write out pickle of the PSSMs for the systematic SELEX experiemnts
outFile = open('selexPSSMs.pkl','wb')
cPickle.dump(selexPSSMs,outFile)
outFile.close()

# 4. Compare motifs to identify non-redundant set
if not os.path.exists('tmp/tomtom'):
    os.makedirs('tmp/tomtom')
outFile = open('upstreamMotifPermutedPValues.csv','w')
outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
pssmsNames = selexPSSMs.keys()
print 'Making files...'
tested = []
for i in range(len(selexPSSMs)):
    print pssmsNames[i]
    if not os.path.exists('tmp/query'+str(i)+'.tomtom'):
        makeQueryFile(nucFreqs={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, queryPssms=[selexPSSMs[pssmsNames[i]]],num=i)
    if not os.path.exists('tmp/target'+str(i)+'.tomtom'):
        targetPSSMs = []
        p1 = pssmsNames[i].split('_')
        if not p1[0].upper()+'_'+p1[3]+'_'+p1[4] in tested:
            tested.append(p1[0].upper()+'_'+p1[3]+'_'+p1[4])
            for pssm in pssmsNames:
                p2 = pssm.split('_')
                if p1[0].upper()==p2[0].upper() and p1[3]==p2[3] and p1[4]==p2[4] and not pssm==pssmsNames[i]:
                    print 'In.'
                    
                    targetPSSMs.append(selexPSSMs[pssm])
        makeTargetFile(nucFreqs={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, targetPssms=targetPSSMs, num=i)
print 'Done.'

# Run this using all cores available
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.' 
#runTomTom(0)
pool = Pool(processes=cpus)
pool.map(runTomTom,range(len(selexPSSMs)))
print 'Done with Tomtom runs.\n'


print 'Removing redundancy in PSSMs...'
redundantPSSMs = []
for run in range(len(selexPSSMs)):
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
            if len(tmpPssms)==0:
                tmpPssms.append(splitUp[0])
            tmpPssms.append(splitUp[1])
        # Decision tree for choosing non-redundant examplar
        # 1. Are there any human TFs?
        mouseRE = re.compile('^[A-Z][a-z]*$')
        human = []
        mouse =  []
        exemplar = []
        for pssm1 in tmpPssms:
            if mouseRE.match(pssm1):
                mouse.append(pssm1)
            else:
                human.append(pssm1)
        # If ther are human motifs
        if len(human)>0:
            # If there is only one human motif
            if len(human)==1:
                redundantPSSMs += mouse
                exemplar = human[0]
            # If there is more than one human motif
            else:
                # Choose full length
                full = []
                dbd = []
                for pssm1 in human:
                    if pssm1.split('_')[2]=='full':
                        full.append(pssm1)
                    else:
                        dbd.append(pssm1)
                # If there is a full length clone motif
                if len(full)>0:
                    # If there is only one
                    if len(full)==1:
                        redundantPSSMs += mouse
                        redundantPSSMs += dbd
                        exemplar = full[0]
                    # If ther are more than one full length clone motif
                    else:
                        print 'Fuck. More than one full length clone motif.'
                else:
                    print 'Fuck. No full length clone motifs.'
        else:
            print 'Fuck. No human motifs.'
        print exemplar, tmpPssms


