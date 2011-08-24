import os
import sys
import csv
import collections

from rSeq.utils.errors import *
from rSeq.utils.misc import Bag

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

    reader  = csv.reader(open(tablePath), delimiter=sep)
    if not headers:
        headers = reader.next()   
    Table   = collections.namedtuple('Table', ', '.join(headers))
    data    = map(Table._make, reader)
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
        
        self._file = open(filePath, 'rU')
        if key:
            self._key = key
        else:
            self._key = lambda x:x[1:].split()[0]
        self.bufferLine = None   # stores next headerLine between records.
        self.joinWith = joinWith
        
    
    def __iter__(self):
        return self
        
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqName,seqStr)"""
        # ++++ Get A Record ++++
        recHead = ''
        recData = []
        # ++++ Check to see if we already have a headerLine ++++
        if self.bufferLine:
            recHead = self.bufferLine
        else:
        # ++++ If not, seek one ++++
            while 1:
                line = self._file.readline()
                if line.startswith('>'):
                    recHead = line
                    break
                elif not line:
                    raise InvalidFileFormatError, "CheckFastaFile: Encountered EOF before any data."
                elif line.strip() == '':
                    continue
                else:
                    raise InvalidFileFormatError, 'CheckFastaFile: The first line containing text does not start with ">".'
        # ++++ Collect recData ++++
        while 1:
            line = self._file.readline()
            if not line:
                raise StopIteration
            elif line.startswith('>'):
                self.bufferLine = line.strip('\n')
                break
            elif not line.startswith('>'):
                recData.append(line.strip('\n'))

        # ++++ Minor Seq Validation ++++
        ## AddHere
        # ++++ Format Rec For Return ++++
        if not recData:
            return None
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
               Please check FastQ file near line #%s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[0])
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file and try again **" % (self._hdSyms[1])
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
    
    def get_next_readSeq(self):
        """Convenience method: calls self.next and returns only the readSeq."""
        try:
            record = self.next()
            return record[1]
        except StopIteration:
            return None
        
    def filter_fastQ_headings(self,filteredPath,key=None):
        """
        Iterates through a fastQ file and writes only those recs
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
        


    
    