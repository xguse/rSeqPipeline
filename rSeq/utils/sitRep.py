#__author__="Gus Dunn"
#__date__ ="$Jan 21, 2011 10:30:16 AM$"

import os
import sys
import datetime

# TODO: option to write all logs to one place and soft link to outDir


class SitRepError(Exception):
    """Generic Error Class for SitRep"""
class ArgumentError(SitRepError):
    """Raised when unexpected values for arguments are found."""

def start_sitrep(outDir, verbose=False):
    """Prep and initilize new stdErr/stdOut streams and open log file."""
    # Inform user of what I do
    print """
    NOTE: The sitRep module has been activated!  This redirects this script's
    stdout and stderr to a central location in order to aid in debugging
    and record keeping.
    
    The full log is located in the following directory:
    
    %s""" % (outDir)
    
    
    # Should I start a new file or append to an existing one?    
    now = datetime.datetime.now().isoformat()
    
    try:
        outDirContents = os.listdir(outDir)
    except OSError:
        os.makedirs(outDir)
        outDirContents = os.listdir(outDir)
    
    outDir = outDir.rstrip('/')
    logFiles = []
    for item in outDirContents:
        if item.endswith(".sitrep.log"):
            logFiles.append(item)
    logFileLen = len(logFiles)
    if logFileLen > 1:
        raise SitRepError("startSitRep() found more than one sitrep log!")
    elif logFileLen == 0:
        outFile = open("%s/%s.sitrep.log" % (outDir, now), "a")
    else:
        outFile = open('%s/%s' % (outDir,logFiles[0]), "a")
    # Inform the file of what command and args were used to call my host script
    callTime = datetime.datetime.now().strftime('%H:%M:%S')
    outFile.write('[NEW CMD] [%s]  [DATE: %s]  %s\n' %
        (callTime,now,' '.join(sys.argv)))
    #if verbose:
        #sys.__stdout__.write('[NEW CMD] [%s]  %s\n' %
        #(callTime,' '.join(sys.argv)))

    # Redirect stdOut/stErr
    sys.stdout = StdOut(outFile,verbose=verbose)
    sys.stderr = StdErr(outFile,verbose=verbose)

class StreamHandler(object):
    """Base class to manage redirected std* streams"""
    def __init__(self,outFile,verbose=False):
        """"""
        if not isinstance(outFile,file):
            raise ArgumentError("sitRep.StreamHandler() requires 'outFile' to be file object!")
        self.out = outFile
        
    def  __getattr__(self,attribute):
        """Provide direct access to the attributes of the file obj
        located at self.out."""
        return self.out.__getattribute__(attribute)
    
    def __del__(self):
        """Clean-up method to make sure that all data is flushed
        to file and it is closed cleanly."""

        try:
            self.out.write('==========\n\n')
            self.flush()
            self.close()
        except ValueError:
            pass

    #def write(self):
        #"""Dummy method to be over-ridden by children"""
    #def flush(self):
        #"""Flush current lines to self.out."""
        #self.out.flush()    
    #def close(self):
        #"""Close self.out."""
        #self.out.close()
        
class StdOut(StreamHandler):
    """Mimic a file obj but include useful info in output."""

    def __init__(self,outFile,verbose=False):
        StreamHandler.__init__(self,outFile,verbose)

    def write(self,text):
        prefix = datetime.datetime.now().strftime('[out]   [%H:%M:%S]')
        for line in text.split('\n'):
            if line:
                self.out.write("%s  %s\n" % (prefix,line))
        self.flush()

class StdErr(StreamHandler):
    """Mimic a file obj but include useful info in output."""

    def __init__(self,outFile,verbose=False):
        StreamHandler.__init__(self,outFile,verbose)

    def write(self,text):
        prefix = datetime.datetime.now().strftime('[ERR]   [%H:%M:%S]')
        for line in text.split('\n'):
            if line:
                self.out.write("%s  %s\n" % (prefix,line))
        self.flush()
