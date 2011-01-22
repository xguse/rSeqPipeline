__author__="Gus Dunn"
__date__ ="$Jan 21, 2011 10:30:16 AM$"

import sys
import datetime


class StreamHandler(object):
    """Base class to manage redirected std* streams"""
    def __init__(self,filePath,runName=None): 
        """"""
        self.startTime = datetime.datetime.now().isoformat()
        self.out = open(filePath,"a")
    
    def __del__(self):
        """Clean-up meathod to make sure that all data is flushed
        to file and it is closed cleanly."""
        self.flush()
        self.close()
        
    def write(self):
        """Dummy method to be over-ridden by children"""
    def flush(self):
        """Flush current lines to self.out."""
        self.out.flush()    
    def close(self):
        """Close self.out."""
        self.out.close()
        
class stdOut(object):
    """Mimic a file obj but include useful info in output."""
    
        

    def write(self,text):
        prefix = datetime.datetime.now().strftime('[%H:%M:%S] ')

class stdErr(object):
    """Mimic a file obj but include useful info in output."""