"""
1: ...
"""

import sys,os
import argparse
import logging
import shutil
from tempfile import NamedTemporaryFile


from rSeq.utils.errors import *
from rSeq.utils.misc import Bag
from rSeq.utils.files import fastaRec_length_indexer,tableFile2namedTuple,mv_file_obj
from rSeq.utils.expression import pearsonExpnFilter,mangle_expn_vectors
from rSeq.utils.externals import mkdirp

def 