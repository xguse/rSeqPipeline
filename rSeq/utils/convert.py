import sys

import pysam

def bam2bed(bamPath,bedFile,region=None):
    """Opens bamPath with pysam to access its contents then converts
    each read's info to bed format and writes to bedFile.

    bamPath: path string
    bedFile: an open file obj or stdOut obj
    region:  samtools region string

    modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
    """
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