import os
import sys
import csv
import collections
import gzip
import shutil
import tempfile


from rSeq.utils.errors import *
from rSeq.utils.misc import Bag,fold_seq
from rSeq.utils.externals import runExternalApp


def mv_file_obj(fileObj,newPath='',chmod=False):
    """
    GIVEN:
    1) fileObj: file object to be moved 
    2) newPath: new path (if not abs path, current working dir is used)
       MUST at LEAST include a name for new file.
    3) chmod: string to use for new file permissions (if False: no change).
    
    DO:
    1) fileObj.flush()
    2) use shutil.move(src,dest) to move file location.
    3) change fileObj.name to the new name.
    4) if file is tmpFile, make sure fileObj.delete == False
       (I assume that if you care enough to move the file, you want to keep it)
    5) chmod on new file to set permissions (remember to use octal: 0755 NOT 755)
    
    RETURN:
    1) fileObj
    """
    if not 'file' in str(fileObj):
        raise SanityCheckError("fileObj does not seem to be a file-like object: %s" % (str(fileObj)))

    
    try:
        fileObj.delete = False
    except AttributeError:
        pass
    if not fileObj.closed:
        fileObj.close()
    
    if newPath == '':
        newPath = fileObj.name.split('/')[-1]
    if (not newPath.startswith('/')) and (not newPath.startswith('./')) and (not newPath.startswith('../')):
        newPath = './%s' % (newPath)
        
    shutil.move(fileObj.name,newPath)
    
    fileObj.name = os.path.abspath(newPath)
    if not chmod == False:
        os.chmod(fileObj.name,int(chmod))
        
    return fileObj
    
def fastaRec_length_indexer(fastaFiles):
    """
    GIVEN:
    1) fastaFiles: list of fasta files or dirs containing fasta files
    
    DO:
    1) iterate through all fasta files recording recName and length to a dict
    
    RETURN:
    1) dict with recName and lengths
    
    NOTES:
    1) will complain if it sees more than one fastaRec with the same name ONLY
       if one of the length values disagrees with those already seen.
    """
    recDict = {}
    seqDict = {}
    tmpDict = collections.defaultdict(list)
    
    for each in fastaFiles:
        try:
            # if each is a directory, measure all recs in all fasta files in that dir (ignore subdirs)
            paths = os.listdir(each)
            for p in paths:
                try:
                    p = ParseFastA(p)
                    for name,seq in p:
                        tmpDict[name].append(len(seq))
                        seqDict[name] = seq
                except IOError:
                    # most likely p was a dir, ignore
                    ## TODO: logging code here to inform when this happens
                    pass
                except InvalidFileFormatError:
                    # most likely p did not have valid fasta format, ignore
                    ## TODO: logging code here to inform when this happens
                    pass
                        
        except OSError as errTxt:
            if not ('Not a directory' in errTxt):
                raise
            else:
                # if each is a file, measure all recs in file
                try:
                    p = ParseFastA(each)
                    for name,seq in p:
                        tmpDict[name].append(len(seq))
                        seqDict[name] = seq
                except InvalidFileFormatError:
                    # most likely p did not have valid fasta format, ignore
                    ## TODO: logging code here to inform when this happens
                    pass
    
    
    for rec,lengths in tmpDict.iteritems():
        if not (len(set(lengths)) == 1):
            # SANITY_CHECK: make sure that any duplicate fastaRecs gave the same length, if not: complain and die
            raise SanityCheckError("Encountered fastaRec with lengths that do not agree: %s:%s" % (rec,lengths))
        else:
            # consolodate the lengths lists to a single number
            recDict[rec] = lengths[0]
    
    return (recDict,seqDict)
    
    

def filter_PEfastQs(filterFunc,fwdMatePath,revMatePath,matchedPassPath1,matchedPassPath2,singlePassPath,nonPassPath):
    """
    Takes the paths to mated PE fastq files with coordinated read-ordering.
    Tests whether paired reads satisfy the provided filterFunc.
    For example fastQs from hudsonAlpha should have either "Y" or "N" flag in their header:
    
    @HWI-ST619:70:B0BMTABXX:3:1102:9652:78621 1:N:0:TAGCTT
    GAATGCGATAGTCACAAGGCATGCCGTTGAATATTCGCAACTGAGCTTCG
    +
    @@@?DAD;DADAD:EGAGH=CFEEHB?GIEBFBGI@CFGD@FG>>@GGGH
    
    The fastQ parsers represent each read as a 4 member list with each line's data as an index.

    
    The filterFunc for this case might be:
    lambda x: x[0].split(' ')[-1].split(':')[1] == "N"
    
    If both mates satisfy the filter, the fwd mate is written to matchedPassPath1 and rev mate to matchedPassPath2.
    If only one mate satisfies the filter, it is written to singlePassPath regardles of fwd/rev.
    All reads that do not satisfy the filter are written to nonPassPath.
    
    Notes:
    * The filterFunc does not have to be a simple lambda, but even something like "testMeanQualScore()",
      as long as it returns a True/False with True meaning that the read should be KEPT.
    * Write-files are overwritten if they exist, created otherwise.
    """
    
    
    fwdMates = ParseFastQ(fwdMatePath)
    revMates = ParseFastQ(revMatePath)
    mPassF_file = open(matchedPassPath1, 'w')
    mPassR_file = open(matchedPassPath2, 'w')
    sPass_file  = open(singlePassPath, 'w')
    nPass_file  = open(nonPassPath, 'w')
    
    outFiles = [mPassF_file,
                mPassR_file,
                sPass_file,
                nPass_file]
    
    counts = Bag({'pairs_passed':0,
                  'fwd_passed_as_single':0,
                  'rev_passed_as_single':0,
                  'fwd_failed':0,
                  'rev_failed':0,
                  'total':0})
    
    while 1:
        # get next fastq Records or set mates to None
        try:
            fwdMate = fwdMates.next()
            counts.total += 1
        except StopIteration:
            fwdMate = None
        
        try:
            revMate = revMates.next()
            counts.total += 1
        except StopIteration:
            revMate = None
        
        # break out if both files are empty 
        if (fwdMate == None) and (revMate == None):
            break

        # set up logical switches to guide the outWriting/counting
        if fwdMate == None:
            keepFwd = None
        else:
            keepFwd = filterFunc(fwdMate)
        if revMate == None:
            keepRev = None
        else:
            keepRev = filterFunc(revMate)
            
        #print "fwdMate %s\nkeepFwd %s\nrevMate %s\nkeepRev %s" % (fwdMate[0],keepFwd,revMate[0],keepRev)
        #break
            
        # write fwd and rev to appropriate files and += respective counts
        if keepFwd and keepRev:
            mPassF_file.write('%s\n' % ('\n'.join(fwdMate)))
            mPassR_file.write('%s\n' % ('\n'.join(revMate)))
            counts.pairs_passed += 1
        if keepFwd and not keepRev:
            sPass_file.write('%s\n' % ('\n'.join(fwdMate)))
            nPass_file.write('%s\n' % ('\n'.join(revMate)))
            counts.fwd_passed_as_single += 1
            counts.rev_failed += 1
        if not keepFwd and keepRev:
            nPass_file.write('%s\n' % ('\n'.join(fwdMate)))
            sPass_file.write('%s\n' % ('\n'.join(revMate)))
            counts.fwd_failed += 1
            counts.rev_passed_as_single += 1
        if not keepFwd and not keepRev:
            nPass_file.write('%s\n' % ('\n'.join(fwdMate)))
            nPass_file.write('%s\n' % ('\n'.join(revMate)))
            counts.fwd_failed += 1
            counts.rev_failed += 1
    
    for f in outFiles:
        f.close()
    
    reportTxt = '''================
Filtered your files using the supplied filter function:
fwdMatePath:\t\t%s
revMatePath:\t\t%s
matchedPassPath1:\t%s
matchedPassPath2:\t%s
singlePassPath:\t\t%s
nonPassPath:\t\t%s\n----\n''' % (fwdMatePath,revMatePath,matchedPassPath1,matchedPassPath2,singlePassPath,nonPassPath)
    sanityCount = (counts.pairs_passed * 2) + counts.fwd_passed_as_single + counts.rev_passed_as_single + counts.fwd_failed + counts.rev_failed
    if not sanityCount == counts.total:
        reportTxt += "WARNING: sanityCount (%s) does not equal counts.total (%s)\n" % (sanityCount, counts.total)
    reportTxt += "PairsPassed:\t\t%s\nFwdSinglePassed:\t%s\nRevSinglePassed:\t%s\nFwdFailed:\t\t%s\nRevFailed:\t\t%s\n" % (counts.pairs_passed,
                                                                                                                           counts.fwd_passed_as_single,
                                                                                                                           counts.rev_passed_as_single,
                                                                                                                           counts.fwd_failed,
                                                                                                                           counts.rev_failed)
    sys.stderr.write("%s\n" % (reportTxt))
    

#def strip_str_of_comments(string,commentStr='#'):
    #"""
    #Take string and return everything to the right of the first
    #comment string found, unless it is quoted or escaped: '#' or "#" or \#.
    #"""
    #qComStr = """'%s'""" % (commentStr)
    #qqComStr = '''"%s"''' % (commentStr)
    #esqComStr = '\%s' % (commentStr)
    #if """'%s'""" % (commentStr) in string:
        
    #return string.rstrip('/n').split(commentStr)[0]

    
def open_wo_comments(filePath,commentStr):
    """
    Return a generator that acts like fileObj.next() but removes
    both commented lines and in-line comments. 
    
    Preserves line numbers.
    """
    inFile = open(filePath,'rU')
    for line in inFile:
        line = strip_str_of_comments(line)
        yield '%s\n' % (line) # replace \n to preserve the behavior of fileObj.next()



def unSoftMask(inFastaPath,outFastaPath):
    # TODO: replace unSoftMask() with mask_converter() that provides multiple options for changing the masking of sequences
    """
    UPPERcases any lowercased nucs in the fasta recs.
    Writes new file.
    """
    inFasta  = open(inFastaPath, 'rU')
    outFasta = open(outFastaPath, 'w')
    for line in inFasta:
        if line.startswith('>'):
            outFasta.write(line)
        else:
            outFasta.write(line.upper())
    inFasta.close()
    outFasta.close()
    
    
def tableFile2namedTuple(tablePath,sep='\t',headers=None):
    """Returns namedTuple from table file using first row fields as
    col headers or a list supplied by user."""

    reader  = csv.reader(open(tablePath,'rU'), delimiter=sep)
    if not headers:
        headers = reader.next()   
    Table   = collections.namedtuple('Table', headers)
    # wrap Table.__getattribute__() for less typing
    def get(self,colName):
        return self.__getattribute__(colName)
    Table.get = get
    
    data    = [Table._make(x) for x in reader if x!=[]] # reader kept feeding an empty list at the end that botched everything!  wtf?!
    return data

class ParseFastA(object):
    """Returns a record-by-record fastA parser analogous to file.readline()."""
    def __init__(self,filePath,joinWith='',key=None):
        """Returns a record-by-record fastA parser analogous to file.readline().
        Exmpl: parser.next()
        Its ALSO an iterator so "for rec in parser" works too!
        
        <joinWith> is string to use to join rec lines with.
        joinWith='' results in a single line with no breaks (usually what you want!)
        
        <key> is func used to parse the recName from HeaderInfo.
        """
        
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
            
        if key:
            self._key = key
        else:
            self._key = lambda x:x[1:].split()[0]
        self.bufferLine = None   # stores next headerLine between records.
        self.joinWith = joinWith
        self._stop = False
        
    
    def __iter__(self):
        return self
        
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqName,seqStr)"""
        if not self._stop:
            pass
        else:
            raise StopIteration
        # ++++ Get A Record ++++
        recHead = ''
        recData = []
        # ++++ Check to see if we already have a headerLine ++++
        if self.bufferLine:
            recHead = self.bufferLine
        else:
        # ++++ If not, seek one ++++
            while 1:
                try:
                    line = self._file.next()
                except StopIteration:
                    self._stop = True
                    break
                if line.startswith('>'):
                    recHead = line
                    break
                elif not line:
                    raise InvalidFileFormatError, "CheckFastaFile: Encountered EOF before any data."
                elif line == '\n':
                    continue
                else:
                    raise InvalidFileFormatError, 'CheckFastaFile: The first line containing text does not start with ">".'
        # ++++ Collect recData ++++
        while 1:
            try:
                line = self._file.next()
            except StopIteration:
                self._stop = True
                break
            if line.startswith('>'):
                self.bufferLine = line.strip('\n')
                break
            elif not line.startswith('>'):
                recData.append(line.strip('\n'))

        # ++++ Minor Seq Validation ++++
        ## AddHere
        # ++++ Format Rec For Return ++++
        if not recData:
            recHead = self._key(recHead)
            return (recHead,'')
        else:
            recHead = self._key(recHead)
            return (recHead,self.joinWith.join(recData))   
    
    def to_dict(self):
        """Returns a single Dict populated with the fastaRecs
        contained in self._file."""
        fasDict = {}
        while 1:
            try:
                fasRec = self.next()
            except StopIteration:
                break
            if fasRec:
                if not fasRec[0] in fasDict:
                    fasDict[fasRec[0]] = fasRec[1]
                else:
                    raise InvalidFileFormatError, "DuplicateFastaRec: %s occurs in your file more than once."
            else:
                break
        return fasDict
    
    def rewrite_headers(self,outPath,lineLen=70,delim=' ',order=[],ow=False,chmod=755):
        """
        PURPOSE
        * reorganize the headers of a fasta file:
          >supercontig:CpipQ1:supercont3.1:1:3873040:1 supercontig supercont3.1
          to 
          >supercont3.1 supercontig:CpipQ1:supercont3.1:1:3873040:1 supercontig
        
        NOTES
        1) If ow, ignores outPath
        2) delim is what to spilt on
        3) order is a list of index numbers from the original header, reorganized for the new header.
           Exp: delim=' ',order=[2,0,1]  would produce what is seen above.
        4) chmod= set new file with this mode (exp: 755)
        5) lineLen sets fastaSeq line length in new file.
        """
        
        if ow:
            outPath = tempfile.NamedTemporaryFile(suffix='.renamed.fas')
        else:
            outPath = open(outPath,'w')
        
        while 1:
            try:
                f = self.next()
            except StopIteration:
                break
            fSplit  = f[0].lstrip('>').rstrip('\n').split(delim)
            newHead = delim.join([fSplit[x] for x in order])
            outPath.write('>%s\n%s\n' % (newHead,'\n'.join(fold_seq(f[1],lineLen))))
        
        self._file.close()
        
        absPath     = os.path.abspath
        if ow:
            outPath.flush()
            outPath.delete = False
            os.rename(absPath(outPath.name),absPath(self._file.name))
                
            try:
                chmodResult = runExternalApp('chmod', '%s %s' % (chmod,absPath(self._file.name)))
            except ExternalError as err:
                sys.stderr.write('%s\n' % (err))
        else:
            try:
                outPath.close()
                chmodResult = runExternalApp('chmod', '%s %s' % (chmod,absPath(outPath.name)))
            except ExternalError as err:
                sys.stderr.write('%s\n' % (err))
        
    



class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
        
    def __iter__(self):
        return self
    
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
        
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
    
    def get_next_readSeq(self):
        """Convenience method: calls self.next and returns only the readSeq."""
        try:
            record = self.next()
            return record[1]
        except StopIteration:
            return None
        
    def filter_SEfastQ_headings(self,filteredPath,key=None):
        """
        Iterates through a single-end fastQ file and writes only those recs
        that satisfy the <key> lambda func to <filteredPath>.
        """
        fastqLen = 0
        filteredLen = 0
        
        if key == None:
            key = lambda x: x
        filtered = open(filteredPath, 'w')
        
        while 1:
            try:
                qRec = self.next()
                fastqLen += 1
            except StopIteration:
                break
            
            if key(qRec[0]) == True:
                for line in qRec:
                    filtered.write('%s\n' % (line))
                    filteredLen += 1
            elif key(qRec[0]) == False:
                pass
            else:
                raise UnexpectedValueError("ERROR: in ParseFastQ.filter_fastQ_headings() 'key' returned a non-T/F value.")
        
        filtered.flush()
        filtered.close()
        
        return Bag({"path":os.path.abspath(filtered.name),
                "original":fastqLen,
                "filtered":filteredLen})
        
def onlyInA(fileA,fileB,outFile):
    """Takes two files. Writes a third file with the lines that are unique to
    fileA.  NOTE: Can be Memory intensive for large files."""
    
    fileA = open(fileA,'rU')
    fileB = open(fileB,'rU')
    
    lineDict = {}
    
    # Get line info
    for line in fileA:
        line = line.strip('\n')
        if line in lineDict:
            lineDict[line][0] += 1
        else:
            lineDict[line] = [1,0]
    
    for line in fileB:
        line = line.strip('\n')
        if line in lineDict:
            lineDict[line][1] += 1
        else:
            lineDict[line] = [0,1]
    
    # write lines unique to fileA
    outFile = open(outFile,'w')
    rLines = 0
    for item in lineDict.iteritems():
        line   = item[0]
        counts = item[1]
        if counts[1] == 0:
            if counts[0] == 1:
                outFile.write('%s\n' % (line))
            elif counts[0] > 1:
                rLines += 1
                outFile.write('%s\n' % (line))
            elif counts[0] == 0:
                raise UnexpectedValueError('It looks like line "%s" does not occur in either file!!\n"lineDict" may have been corrupted.' % \
                                           (line))
        else:
            continue    
    if rLines:
        sys.stderr.write('WARNING: %s lines written to %s occured more than once in %s, but %s is now non-redundant.' %\
                         (rLines,outFile.name,fileA.name,outFile.name))
    outFile.close()
        


    
def renameChrom_in_SAM(path):
    """
    """
    import fileinput
    path = os.path.abspath(path)
    bkExt = ".zap_me.backup"
    bkPath = path + bkExt
    
    try:
        for line in fileinput.input(files=path, inplace=1, backup=bkExt,
                                    bufsize=0,mode="r", openhook=None):
            if line.startswith('@SQ'):
                line = line.split('\t')
                chrm = "SN:%s" % (line[1].split(':')[3])
                line[1] = chrm
                sys.stdout.write('\t'.join(line))
            elif line.startswith('HWI'):
                line = line.split('\t')
                chrm = line[2].split(':')[2]
                line[2] = chrm
                sys.stdout.write('\t'.join(line))
            else:
                sys.stdout.write(line)
    except:
        raise


def rename_fasta_headers(in_path,out_path,header_func):
    """
    GIVEN:
        - in_path = path to original fasta file
        - out_path = path to future altered fasta file
        - header_func = function to take a header line and return an altered string version of it
    DOES:
        - Reads in in_path file one line at a time
        - If the line is a fasta header (starts with '>')
          uses header_func logic to rearrange the header and
          writes out the changed line to out_path.
        - If not a header, writes same line out to out_path.
        - Closes both file objects.
    RETURNS:
        - None
    """
    
    in_file = open(in_path,'rU')
    out_file = open(out_path,'w')
    
    for line in in_file:
        if line.startswith('>'):
            line = header_func(line)
            # Handle and ensure that each modified line has one and only one \n
            line = line.rstrip('\n') + '\n' 
        else:
            pass
        
        out_file.write(line)
    
    in_file.close()
    out_file.close()
