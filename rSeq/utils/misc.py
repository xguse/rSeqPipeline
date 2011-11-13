import sys
import textwrap
from optparse import IndentedHelpFormatter


class RseqHelpFormatter(IndentedHelpFormatter):
    """Custom help formatter for optparse that allows
    me to break text into paragraphs!"""
    
    def __init__(self,
                 indent_increment=0,
                 max_help_position=80,
                 width=None,
                 short_first=1):
        IndentedHelpFormatter.__init__(self,
                                       indent_increment=indent_increment,
                                       max_help_position=max_help_position,
                                       width=width,
                                       short_first=short_first)
    def _format_text(self, text):
        """
        Format a paragraph of free-form text for inclusion in the
        help output at the current indentation level while honoring
        users' newlines!
        """
        text_width = self.width - self.current_indent
        indent = " "*self.current_indent
        return textwrap.fill(text,
                             text_width,
                             initial_indent='',
                             subsequent_indent='\t',
                             replace_whitespace=False,
                             drop_whitespace=True,
                             expand_tabs=True)
    


# From Titus Brown's gff parser:
class Bag(dict):
    """dict-like class that supports attribute access as well as getitem.

    >>> x = Bag()
    >>> x['foo'] = 'bar'
    >>> x.foo
    'bar'

    """
    def __init__(self, *args, **kw):
        dict.__init__(self, *args, **kw)
        for k in self.keys():
            self.__dict__[k] = self.__getitem__(k)

    def __setitem__(self, k, v):
        dict.__setitem__(self, k, v)
        self.__dict__[k] = v
        
def pVal4mari(tabPath,tTx,tTxDN,tTxUP):
    """
    tTx=   total transcripts,
    tTxUP= total Tx sig up Reg
    tTxDN= total Tx sig dwn reg,
    
    Takes named tuple data table in format:
    type,total,down,up,tDE
    
    Returns list of pValues in format:
    [(type1,pVal_dwn1,pVal_up1,pVal_tde1),
    ...,
    (typeN,pVal_dwnN,pVal_upN,pVal_tdeN)]
    """
    from scipherSrc.defs.files_io import tableFile2namedTuple
    from scipy.stats import hypergeom
    def cHgPvl(x,M,n,N):
        """
        x=randVar
        M=popSize
        n=totalSuccesses
        N=samplSize
        """
        return 1-hypergeom.cdf(x,M,n,N)+hypergeom.pmf(x,M,n,N)
    
    namdTup = tableFile2namedTuple(tabPath,sep=',')
    rData = []
    for i in range(len(namdTup)):
        t = namdTup[i]
        pDwn = cHgPvl(int(t.down),tTx,int(t.total),tTxDN)
        pUp  = cHgPvl(int(t.up),tTx,int(t.total),tTxUP)
        pDE  = cHgPvl(int(t.tDE),tTx,int(t.total),(tTxDN+tTxUP))
        rData.append((t.type,pDwn,pUp,pDE))
    return rData
        
        
    
def whoami():
    """Returns the name of the currently active function."""
    return sys._getframe(1).f_code.co_name

def slidingWindow(sequence,winSize,step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
    
    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
    
    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize)/step)+1
    
    # Do the work
    for i in range(0,numOfChunks,step):
        yield sequence[i:i+winSize]