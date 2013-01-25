import cPickle
from pssm import pssm
from copy import deepcopy

# Make the files for a TomTom run
def makeFiles(nucFreqs, queryPssms, targetPssms, num, strands='+ -'):
    # Header crap
    memeHeader = ''
    memeHeader += 'MEME version 3.0\n\n'
    memeHeader += 'ALPHABET= ACGT\n\n'
    # Here is where we tell it what strand: for miRNAs this would just be '+'
    memeHeader += 'strands: '+strands+'\n\n'
    memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
    memeHeader += 'A '+str(round(float(nucFreqs['A']),3))+' C '+str(round(float(nucFreqs['C']),3))+' G '+str(round(float(nucFreqs['G']),3))+' T '+str(round(float(nucFreqs['T']),3))
    # Make query PSSM file
    queryFile = open('tmp/query'+str(num)+'.tomtom','w')
    queryFile.write(memeHeader)
    queryFile.write('\n\n'.join([pssm1.getMemeFormatted() for pssm1 in queryPssms]))
    queryFile.close()
    # Make target PSSM file
    targetFile = open('tmp/target'+str(num)+'.tomtom','w')
    targetFile.write(memeHeader)
    targetFile.write('\n\n'.join([pssm1.getMemeFormatted() for pssm1 in targetPssms]))
    targetFile.close()

# Run TomTom on the files
def TomTom(num, distMeth='ed', qThresh='1', minOverlap=6):
    # Arguments for tomtom
    tomtomArgs = ' -query tmp/query'+str(num)+'.tomtom -target tmp/target'+str(num)+'.tomtom -dist '+str(distMeth)+' -o tmp/tomtom_out -text -q-thresh '+str(qThresh)+' -min-overlap '+str(minOverlap)+' -verbosity 0'
    print tomtomArgs
    #p = Popen("tomtom" + tomtomArgs, shell=True)
    #sts = os.waitpid(p.pid, 0)
    errOut = open('tmp/stderr.out','w')
    tomtomProc = Popen("tomtom" + tomtomArgs, shell=True,stdout=PIPE, stderr=errOut)
    outputFile = open('tmp/tomtom_out/tomtom'+str(num)+'.out','w')
    output = tomtomProc.communicate()[0]
    outputFile.write(output)
    outputFile.close()
    errOut.close()

# Wrapper function to run TomTom using multiprocessing pool
def runTomTom(i):
    TomTom(i, distMeth='ed', qThresh='1', minOverlap=6) #blic5

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
outFile = open('upstreamMotifPermutedPValues.csv','w')
outFile.write('Motif Name,Region,Original E-Value,Consensus,Permuted E-Value < 10,Similar,Total Permutations,Permuted P-Value')
pssmsNames = pssms.keys()
print 'Making files...'
for i in range(len(selexPSSMs)):
    if not os.path.exists('tmp/query'+str(i)+'.tomtom') and not os.path.exists('tmp/target'+str(i)+'.tomtom'):
        makeFiles(nucFreqs={'A':0.25,'C':0.25,'G':0.25,'T':0.25}, queryPssms=[pssms[pssmsNames[i]]],targetPssms=selexPSSMs,num=100)
print 'Done.'

# Run this using all cores available
cpus = cpu_count()
print 'There are', cpus,'CPUs avialable.' 
pool = Pool(processes=cpus)
pool.map(runTomtom,range(len(pssms)))
print 'Done with Tomtom runs.\n'

print 'Reading in Tomtom run...'
for run in range(len(pssms)):
    tomtomPValues = {}
    outputFile = open('tmp/tomtom_out/tomtom'+str(run)+'.out','r')
    output = outputFile.readlines()
    outputFile.close()
    # Now iterate through output and save data
    output.pop(0) # Get rid of header
    while len(output)>0:
        outputLine = output.pop(0).strip().split('\t')
        if len(outputLine)==9:
            tomtomPValues[outputLine[1]] = float(outputLine[3])
    pValues = tomtomPValues.values()
    similar = 0
    for pValue in pValues:
        if float(pValue) <= float(0.05):
            similar += 1
    # Write out the results
    mot = outputLine[0].split('_')[1]
    permPValues[outputLine[0]] = { mot+'.consensus':str(pssms[outputLine[0]].getConsensusMotif()), mot+'.permutedEV<=10':str(len(pValues)), mot+'.similar':str(similar), mot+'.permPV':str(float(similar)/float(1000)) }
    if outputLine[0] in upstreamMatches.keys():
        permPValues[outputLine[0]]['motif1.matches'] = ' '.join(upstreamMatches[outputLine[0]])
        matched += 1
    else:
        permPValues[outputLine[0]][mot+'.matches'] = 'NA'
    outFile.write('\n'+str(outputLine[0])+',upstream,'+str(pssms[outputLine[0]].getEValue())+','+str(pssms[outputLine[0]].getConsensusMotif())+','+str(len(pValues))+','+str(similar)+','+str(1000)+','+str(float(similar)/float(1000)))
outFile.close()

