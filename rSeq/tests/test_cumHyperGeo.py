from rSeq.utils.stats import cumHypergeoP

    #Calculates the cumulative hypergeometric p-value for variables:
    #n = # of positives in population
    #i = # of positives in sample
    #m = # of negatives in population
    #N = sample size

    #P(i) = sum([as i->N] (choose(n,i)choose(m,N-i))/choose(n+m,N))

    #For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html


None