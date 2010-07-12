import os
import random
from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp,runExternalApp

def runSCOPE(pLen,genes,jobName,scopeDir,outDir,paramName,jMem='2000',verbose=False):
    """Perform a SCOPE run. Complain and quit if error occurs.
    Notes:
    scopeDir = full path.
    genes    = 'gene;gene;gene;etc'
    pLen     = promorter length to use."""
    # Get full path (if not given) for outDir since we will be jumping around in the directory tree
    if not outDir.startswith('/'):
        outDir = os.getcwd()+outDir
        outDir = outDir.rstrip('/')
    else:
        outDir = outDir.rstrip('/')
    
    # Set up argString
    outPathBase = '%s/%s.%s' % (outDir,jobName,pLen)
    argString = '''-Xmx%sm -cp dist/scope.jar edu.dartmouth.bglab.beam.CGIScope -pf "%s" -ofx "%s.xml" -oft "%s.txt" -oje "%s" -qg "%s" -sgl "%s" -drb "true" -dra "true" -drbp "true"''' \
              % (jMem,
                 paramName,
                 outPathBase,
                 outPathBase,
                 jobName,genes,pLen)
    
    # Change to scopeDir for execution bc SCOPE is a PITA.
    os.chdir(scopeDir)
    mkdirp(outDir) # make outDir along with parent dirs as needed
    print 'starting run...'
    resultSCOPE = runExternalApp('java',argString)
    
    # write stdOut/Err to files if requested
    if verbose:
        stdOutFile = open(outPathBase+'.out','w')
        stdErrFile = open(outPathBase+'.err','w')
        stdOutFile.write(resultSCOPE[0])
        stdErrFile.write(resultSCOPE[1])
        stdOutFile.close()
        stdErrFile.close()
    
    return resultSCOPE
        
def listParams(scopeHomeDir):
    """Print list of param files in given scope HomeDir."""
    scopeHomeDir = scopeHomeDir.rstrip('/')
    for x in os.listdir(scopeHomeDir+'/data/params/'):
        print x
        
def getRandomGeneSet(scopeHomeDir,paramName,number):
    """Return list of randomly sampled genes of length number cooresponding
    to paramSet: paramName."""
    scopeHomeDir = scopeHomeDir.rstrip('/')
    paramName = paramName.rstrip('.param')
    genePop = map(lambda l: l.split()[0], open('%s/data/geneLists/%s_IDs.txt' % (scopeHomeDir,paramName), 'rU'))
    
    return random.sample(genePop,number)
    