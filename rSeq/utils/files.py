import sys
import csv
import collections
from rSeq.utils.errors import *

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
    
    def toDict(self):
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
        
        