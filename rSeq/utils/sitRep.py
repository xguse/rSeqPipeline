__author__="Gus Dunn"
__date__ ="$Jan 21, 2011 10:30:16 AM$"

import sys
import datetime


class stdOut(object):
    """Mimic a file obj but include useful info in output."""
    def __init__(self,fileObj): # file should be appendable
        self.out = fileObj

    def write(self,text):
        prefix = datetime.datetime.now().strftime('[%H:%M:%S] ')

class stdErr(object):
    """Mimic a file obj but include useful info in output."""