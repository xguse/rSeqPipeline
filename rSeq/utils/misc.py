import sys

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
