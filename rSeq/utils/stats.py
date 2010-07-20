from decimal import Decimal
# see bottom for conditional import of "bestChoose"



def benjHochFDR(pVals,pValColumn=1,FDR=0.05):
    """
    pVals      = 2D list(hypothesis,p-value) hypothesis could = geneName tested for enrichment
    pValColumn = integer of column index containing the p-value.
    !*! FDR        = threshold above which FDR is unacceptable [not yet implemented]!*! 
    
    Returns, for all *acceptable q-values: hypothesis,origPval,adjustedPval 
    *NOTE:  For now, returns _ALL_ items passed to it with no filtering at the moment.
    """
    assert type(pValColumn) == type(1),\
           "ERROR: pValColumn must be int type!"
    # Sort pVals from highest to lowest after converting them all to floats.
    for i in range(len(pVals)):
        pVals[i][pValColumn] = float(pVals[i][pValColumn])
    pVals.sort(key=lambda x: x[pValColumn])
    pVals.reverse()
    
    n = len(pVals)
    
    lastPval = pVals[0][pValColumn]
    for i in range(len(pVals)):
        p    = pVals[i][pValColumn]
        adj  = (float(n)/(n-i))
        adjP = p*adj
        miN  = min(adjP,lastPval)
        pVals[i].append(miN)
        lastPval = pVals[i][-1]
    
    pVals.reverse()
    return pVals




def binComb(n, k):
    """
    binComb(n, k): Computes n choose k. Defines the number of k objects that can be chosen from 
    n objects.
    """
    if (k > n): return 0
    if (k < 0): return 0
    if (k > int(n/2)):
        k = n - k
    rv = 1
    for j in range(0, k):
        rv *= n - j
        rv /= j + 1
    return rv


def choose(n, k):
    if (k > n): return 0
    if (k < 0): return 0
    ntok = 1
    for t in xrange(min(k, n-k)):
        ntok = ntok*(n-t)//(t+1)
    return ntok

# determine whether we can use gmpy and set "bestChoose" to that or the standby.
try:
    from gmpy import comb as bestChoose #timeit says gmpy.comb is 10 times faster than choose!
except ImportError:
    bestChoose = choose # timeit says that choose is 7% faster than binComb
    
    

def hypergeoP(n,i,m,N):
    """
    Calculates the non-cumulative hypergeometric p-value for variables:
    n = # of positives in population
    i = # of positives in sample
    m = # of negatives in population
    N = sample size

    P(x=i) = (choose(n,i)choose(m,N-i))/choose(n+m,N)

    For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    """
    return (bestChoose(n,i)*bestChoose(m,N-i))/float(bestChoose(n+m,N))




def cumHypergeoP(n,i,m,N):
    """
    Calculates the cumulative hypergeometric p-value for variables:
    n = # of positives in population
    i = # of positives in sample
    m = # of negatives in population
    N = sample size

    P(i) = sum([as i->N] (choose(n,i)choose(m,N-i))/choose(n+m,N))

    For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    """

    cumPVal = 0

    for x in range(i,N+1):
        cumPVal = cumPVal + hypergeoP(n,x,m,N)

    return cumPVal


def binomialPval(n,k,p):
    """Returns exact binomial P-value.
    n = number of trials
    k = number of successes
    p = probability of a success
    
    P(k succeses in n trials) = choose(n,k) * p^k * (1-p)^(n-k)
    """
    
    return bestChoose(n,k) * p**k * (1-p)**(n-k)

def binomialPval_gte(n,k,p):
    """Returns binomial P-value of k or greater successes in n trials."""
    cumPVal = 0
    for k in range(k,n+1):
        cumPVal = cumPVal + binomialPval(n,k,p)
    return cumPVal




def seqStats(seqDict):
    """Returns Dict of useful sequence statistics."""
    combinedSeq = ''
    
    for each in seqDict:
        combinedSeq += seqDict[each]
    
    combinedSeq= combinedSeq.upper()
    
    seqs       = len(seqDict)
    totNucs    = len(combinedSeq)
    aCnt       = combinedSeq.count('A')
    cCnt       = combinedSeq.count('C')
    gCnt       = combinedSeq.count('G')
    tCnt       = combinedSeq.count('T')
    nCnt       = combinedSeq.count('N')
    nonNs      = aCnt+cCnt+gCnt+tCnt
    n2tot      = float(nCnt)/len(combinedSeq)
    n2nonN     = float(nCnt)/nonNs
    percentGC  = (float(gCnt)+cCnt)/nonNs
    
    
    
    return {'seqNum':seqs,
            'totNucs':totNucs,
            'aCnt':aCnt,
            'cCnt':cCnt,
            'gCnt':gCnt,
            'tCnt':tCnt,
            'nCnt':nCnt,
            'nonNs':nonNs,
            'n2tot':n2tot,
            'n2nonN':n2nonN,
            'percentGC':percentGC}