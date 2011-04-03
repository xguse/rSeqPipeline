import sys, pdb

import pysam

from rSeq.utils.errors import *

# import pp if avail
try:
    import pp
except ImportError:
    pass




def bam2bed(bamPath,bedFile,region=None,ncpus="autodetect"):
    """Opens bamPath with pysam to access its contents then converts
    each read's info to bed format and writes to bedFile.

    bamPath: path string or list of path strings
    bedFile: an open file obj or stdOut obj or a list of open file objs
    region:  samtools region string or list of region strings
    ncpus: num of processors to use if parallel job is run (autodetect -> all avail)
    
    If list is given for bamPath only, all BAMs are written to the same
    bedFile
    
    If list given for bamPath AND bedFile, they must be of equal length,
    and bamPaths[i] is written to bedFile[i].
    
    If list is given for region, it must be of same length as bamPath and
    bamPaths[i] will use regions[i].
    
    When appropriate, the jobs will be run in parallel.

    modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
    """
    pdb.set_trace()
    def convertBAM(bamPath,bedFile,region=None):
        # Check to see if the file is closed if so, open it for appending
        #pdb.set_trace()
        if bedFile.closed:
            bedFile = open(bedFile.name,'a')
        # Report convertion files:
        print "Converting %s to %s..." % (bamPath,bedFile.name) 
        # open BAM and get the iterator
        samfile = pysam.Samfile(bamPath,"rb")
        if region != None:
            it = samfile.fetch(region=region)
        else:
            it = samfile.fetch()
    
        # calculate the end position and print out BED
        take = (0, 2, 3) # CIGAR operation (M/match, D/del, N/ref_skip)
        outfile = bedFile
        for read in it:
            if read.is_unmapped: continue
            # compute total length on reference
            t = sum([l for op,l in read.cigar if op in take])
            if read.is_reverse: strand = "-"
            else: strand = "+"
            outfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %\
                          (samfile.getrname(read.rname),
                           read.pos, read.pos+t, read.qname,
                           read.mapq, strand))
            if outfile.name != "<stdout>":
                outfile.flush()
            else:
                outfile.close()
           
    # # Do some validation and set up # #
    
    # convert all that are allowed to be lists into lists
    if type(bamPath) == type(''):
        bamPath = [bamPath]
    if type(bedFile) == "<type 'file'>":
        bedFile = [bedFile]
    if region:
        if region == type(''):
            region = [region]
    else:
        region = [None]
        
    
    # enforce len(bamPathList) == len(bedFileList) unless len(bedFileList) == 1
    if (len(bedFile) != 1) and (len(bamPath) != len(bedFile)):
        raise ValueError('ERROR: If list given for bamPath AND bedFile, they must be of equal length,and bamPaths[i] is written to bedFile[i]. ')
                
    # enforce len(bamPathList) == len(regionList) unless len(regionList) == 1
    if (len(region) != 1) and (len(bamPath) != len(region)):
            raise ValueError('ERROR: If list is given for region, it must be of same length as bamPath andbamPaths[i] will use regions[i]. ')
    
    # enforce type(bedFile) == "<type 'file'> if only one bedFile given
    if (len(bedFile) != 1) and (len(bamPath) != len(bedFile)):
        raise TypeError('ERROR: If list given for bamPath AND NOT bedFile, bedFile must be an open file obj. ')
        
        
    # enforce type(bedFile) == type([]) if only one bamPath given
    if (len(bamPath) == 1) and (len(bedFile) > 1):
        raise TypeError('ERROR: If single bamPath, bedFile must be a single open file obj. ')
    
        
    # equalize list lengths
    if len(bedFile) == 1:
        bedFile = bedFile * len(bamPath)
    if len(region) == 1:
        region = region * len(bamPath)
    
    # # Initialize job_server if pp is avail and run jobs # #
    try:
        # if pp: run in parallel
        job_server = pp.Server(ncpus=ncpus)
        print "Running with %s cpus." % (job_server.get_ncpus())
        jobs = []
        results = []
	#pdb.set_trace()
        for i in range(len(bamPath)):
            jobs.append(job_server.submit(func=convertBAM,
                                          args=(bamPath[i],open(bedFile[i],'a'),region[i]),
                                          depfuncs=(),
                                          modules=('pysam',),
                                          callback=None,
                                          callbackargs=(),
                                          group='default',
                                          globals=None))
	
        for i in range(len(bamPath)):
            results.append(jobs[i]())

        for i in range(len(bamPath)):
            print results[i]

    except NameError as err:
        # if no pp: run sequencially
        if err[0] == "name 'pp' is not defined":
            print "Alert: pp module not found, not running parallel."
            for i in range(len(bamPath)):
                convertBAM(bamPath[i],bedFile[i],region[i])
        


def genephred2refflat(genePhredPath,refFlatPath):
    genePhred = open(genePhredPath,'rU')
    refFlat   = open(refFlatPath,'w')
    for line in genePhred:
        fields = line.rstrip('\n').split('\t')
        for i in range(len(fields)):
            fields[i] = fields[i].rstrip(',')
        if fields[0][-3] == "-":
            fields.insert(0,fields[0][:-3])
        else:
            fields.insert(0,fields[0])
        line = '\t'.join(fields)
        refFlat.write('%s\n' % (line))
    genePhred.close()
    refFlat.close()
