from rSeq.utils.errors import *
from rSeq.utils.externals import runExternalApp,mkdirp
from rSeq.utils.files import onlyInA
def isBEDline(line):
    """Returns TRUE if line 'looks' like a MINIMAL BED formated line."""
    line = line.strip('\n').split('\t')
    if not len(line) >= 3:
        return False
    
    chrm,chrmStart,chrmEnd = line[0],line[1],line[2]
    if chrm.isdigit():
        return False
    elif not chrmStart.isdigit():
        return False
    elif not chrmEnd.isdigit():
        return False
    else:
        return True

def isStranded(line):
    """Returns TRUE if BED line has '+' or '-' at line[5]."""
    line = line.strip('\n').split('\t')
    if not len(line) >= 6:
        return False
    elif not line[5] in ['+','-']:
        return False
    else:
        return True

def extractFromDoubleSidedBedtoolOut(filePath,cols,side='right',outDir='.'):
    """Creates new file from filePath using only the bedInfo from the
    left/right (based on 'side') side of a BEDtools outFile with double-
    sided output (side=[3,6]). 'cols' must be a list with length of columns
    in each 'side' of the double output.  'side' = keep the 'right' or
    'left' side of the output line."""
    
    # Prepare outDir if it doesnt already exist
    mkdirp(outDir)
    
    inFile      = open(filePath, 'rU')
    outFilePath = '%s/%s_%s.bed' % (outDir,
                                    filePath.replace('.bed','').split('/')[-1],
                                    side)
    outFile     = open(outFilePath , 'w')
    lineNum = 0
    for line in inFile:
        lineNum+=1
        if line.startswith('track'):
            continue
        line = line.strip('\n').split('\t')
        
        # Divide the line into two based on cols
        divLine = line[:cols[0]],line[-(len(line)-cols[0]):]

        
        # Ensure the length of each new line is what we expect, then write out cleaned line
        if not ((len(divLine[0])==cols[0]) and (len(divLine[1])==cols[1])):
            raise InvalidFileFormatError('ERROR: line %s in file %s has unexpected number of columns or the values in "cols" is incorrect.' % \
                                         (lineNum,filePath))
        else:
            if side == 'right':
                outFile.write('%s\n' % ('\t'.join(divLine[1])))
            elif side == 'left':
                outFile.write('%s\n' % ('\t'.join(divLine[0])))
            else:
                raise InvalidOptionError('ERROR: option "side" must be one of %s. Was: %s.' % \
                                         (['right','left'], side))
    
    outFile.close()
    # Sort and  remove redundancy from line in new file
    resultSort = runExternalApp('sort','-u %s > %s.tmp' % (outFilePath,outFilePath))
    resultMv   = runExternalApp('mv','%s.tmp %s' % (outFilePath,outFilePath))
    return outFilePath
    
    
def divByWindow(bedA_Path,bedB_Path,win=[500,500],cols=[6,6],side='right',outDir='.'):
    """Create files separating features in bedB by those alling within the area defined
    by <win> and those outside this area in bedA.  If A.bed is stranded, the area is defined
    by win[0] upstrm and win[1] dwnstrm on the FEATURE's strand.  Otherwise its 
    win[0] upstrm and win[1] dwnstrm on the CONTIG/CHROM's plus strand.  Files ouput
    to outDir.
    
    NOTE: See DOC for extractFromDoubleSidedBedtoolOut() regarding 'cols' and 'side'"""

    # Prepare outDir if it doesnt already exist
    mkdirp(outDir)
    
    # Collect some useful info
    bedA_name           = bedA_Path.split('/')[-1].replace('.bed','')
    bedB_name           = bedB_Path.split('/')[-1].replace('.bed','')
    B_in_A_winComboPath = '%s/%s_featsIn_%s_Win%sl%sr_combo.bed' % (outDir,
                                                                    bedB_name,
                                                                    bedA_name,
                                                                    win[0],
                                                                    win[1])
    
    
    # Establish whether inputs look like BED files:
    testA = open(bedA_Path,'rU')
    testB = open(bedB_Path,'rU')
    linesA = []
    linesB = []
    for i in range(2):
        linesA.append(testA.readline())
        linesB.append(testB.readline())
    testA.close()
    testB.close()
    
    if not isBEDline(linesA[1]):
        raise InvalidFileFormatError('ERROR: %s does not seem to be in BED format.' % (bedA_Path))
    if not isBEDline(linesB[1]):
        raise InvalidFileFormatError('ERROR: %s does not seem to be in BED format.' % (bedB_Path))
    
    # If bedA is stranded: use windowBed with -sw option, otherwise with only -l,-r options
    # to create file from bedB features INSIDE window around features in bedA.

    if isStranded(linesA[1]):
        resultWinBed = runExternalApp('windowBed',
                                      '-a %s -b %s -l %s -r %s -sw > %s' % \
                                      (bedA_Path,
                                       bedB_Path,
                                       win[0],
                                       win[1],
                                       B_in_A_winComboPath))
    else:
        resultWinBed = runExternalApp('windowBed',
                                      '-a %s -b %s -l %s -r %s > %s' % \
                                      (bedA_Path,
                                       bedB_Path,
                                       win[0],
                                       win[1],
                                       B_in_A_winComboPath))
    
    # Clean B_in_A_winComboPath of the matching bedA entry and remove any redundant bedB entries
    cleanedBsInWinPath = extractFromDoubleSidedBedtoolOut(B_in_A_winComboPath,
                                                          cols=cols,
                                                          side=side,
                                                          outDir=outDir)
    # Change file name to reflect its not combo anymore
    cleanedBsInWinNewPath = cleanedBsInWinPath.replace('_combo_','_cleaned_')
    resultMv = runExternalApp('mv','%s %s' % \
                              (cleanedBsInWinPath,
                               cleanedBsInWinNewPath))
    
    
    
    # Create file with bedB feats OUTSIDE of window of features in bedA.
    cleanedBsNotInWinPath = cleanedBsInWinNewPath.replace('_featsIn_','_featsNotIn_')
    onlyInA(bedB_Path,cleanedBsInWinNewPath,'%s/%s' % (outDir,cleanedBsNotInWinPath))
    #resultIsectBed = runExternalApp('intersectBed','-a %s -b %s -v > %s' % \
                                    #(bedB_Path,
                                     #cleanedBsInWinNewPath,
                                     #cleanedBsNotInWinPath))
    
    # Return Filenames of divided bed files
    return (cleanedBsInWinNewPath,cleanedBsNotInWinPath)
    