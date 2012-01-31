import sys, pdb

import pysam

from rSeq.utils.errors import *
from rSeq.utils.files import tableFile2namedTuple
from rSeq.utils.align import parseCigarString

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
    #pdb.set_trace()
    def convertBAM(bamPath,bedFile,region=None):
        # Check to see if the file is closed if so, open it for appending
        #pdb.set_trace()
        bedFile = open(bedFile,'a')
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
    if type(bedFile) == type(""):
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

    ## enforce type(bedFile) == "<type 'file'> if only one bedFile given
    #if (len(bedFile) != 1) and (len(bamPath) != len(bedFile)):
        #raise TypeError('ERROR: If list given for bamPath AND NOT bedFile, bedFile must be an open file obj. ')


    # enforce type(bedFile) == type([]) if only one bamPath given
    if (len(bamPath) == 1) and (len(bedFile) > 1):
        raise TypeError('ERROR: If single bamPath, bedFile must be a single bed path. ')


    # equalize list lengths
    if len(bedFile) == 1:
        bedFile = bedFile * len(bamPath)
    if len(region) == 1:
        region = region * len(bamPath)

    ## open files
    #for i in range(len(bamPath)):
    #    bedFile[i] = open(bedFile[i],'a')
        
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
                                          args=(bamPath[i],bedFile[i],region[i]),
                                          depfuncs=(),
                                          modules=('pysam','pdb'),
                                          callback=None,
                                          callbackargs=(),
                                          group='default',
                                          globals=None))

        for i in range(len(bamPath)):
            results.append(jobs[i]())


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
    
def MB_2_gff3(resultTablePath,gff3Path):
    """
    """
    gff3_lines = []
        
    mb_table = tableFile2namedTuple(resultTablePath,sep='\t')
    
    for line in mb_table:
        gff3_seqid = align_feat.chr
        gff3_source = 'Cufflinks'
        gff3_type = 'Assembled Tx boundries'
        gff3_start = left
        gff3_end = right
        gff3_score = line.q_value
        gff3_strand = strandConvertions[align_feat.seq_region_strand]
        gff3_phase = '.'
        gff3_attributes = 'ID=%s;Alias=%s' % (align_feat.dna_align_feature_id, align_feat.hit_name)
        
        gff3_lines.append([gff3_seqid,
                           gff3_source,
                           gff3_type,
                           gff3_start,
                           gff3_end,
                           gff3_score,
                           gff3_strand,
                           gff3_phase,
                           gff3_attributes])        
    

def vectorBaseESTs_2_gff3(resultTablePath,gff3Path):
    """
    GIVEN:
    
    1) resultTablePath: table resulting from mySQL query of EBI's "other features"
       database schema with following query:
           SELECT s.name AS "chr", d.* FROM dna_align_feature d
           LEFT JOIN seq_region s ON (d.seq_region_id = s.seq_region_id)
           WHERE s.seq_region_id < 8 AND d.analysis_id = 227 AND
           s.seq_region_id = d.seq_region_id
    
    2) gffPath: path to new gff3 file
    
    DO:
    
    1) For each row, create multiline GFF3 entries using dna_align_feature_id
       as 'ID' and hit_name as 'Alias', score as 'Score' column, etc. using the 
       cigar_line to generate the feature coords.
    2) Write data to gffPath.
    
    RETURN:
    
    1) Full path of gffPath.
    """
    # helper defs
    def match(current_loc,value):
        """
        """      
        left  = current_loc + 1
        right = left + value - 1
        
        current_loc = right 
        
        return (current_loc,left,right)
    
    def insertion(current_loc,value):
        """
        """
        left  = current_loc + 1
        right = left + value - 1
        current_loc = right 
        
        return (current_loc, None, None)
        
    def deletion(current_loc,value):
        """
        """
        current_loc = current_loc
        
        return (current_loc, None, None)
    
    cigar_operations = {'I':insertion,
                        'M':match,
                        'D':deletion}
    strandConvertions = {'1':'+',
                         '-1':'-'}
    
    gff3_lines = []
    
    est_table = tableFile2namedTuple(resultTablePath,sep=',')
    
    for align_feat in est_table:
        cig_tupl = parseCigarString(align_feat.cigar_line,kind='EBI')
        #cig_tupl = parseCigarString("5M3D6M2I3M",kind='EBI')
        align_feat_lines = []
        
        far_left  = int(align_feat.seq_region_start)
        far_right = int(align_feat.seq_region_end)
        #far_left  = 1
        #far_right = 16        
        
        current_loc = far_left - 1 
        
        for op,value in cig_tupl:
            current_loc, left, right = cigar_operations[op.upper()](current_loc,int(value))
            if left != None:
                #pass
                # Construct the gff line for the match feature and append it to gff_lines
                gff3_seqid = align_feat.chr
                gff3_source = 'Exonerate'
                gff3_type = 'EST_match'
                gff3_start = left
                gff3_end = right
                gff3_score = align_feat.score
                gff3_strand = strandConvertions[align_feat.seq_region_strand]
                gff3_phase = '.'
                gff3_attributes = 'ID=%s;Alias=%s' % (align_feat.dna_align_feature_id, align_feat.hit_name)
                
                align_feat_lines.append([gff3_seqid,
                                         gff3_source,
                                         gff3_type,
                                         gff3_start,
                                         gff3_end,
                                         gff3_score,
                                         gff3_strand,
                                         gff3_phase,
                                         gff3_attributes])
            else:
                pass
            
        # sanity check after each align_feat: does the final location == far_right?
        if not current_loc == far_right:
            raise SanityCheckError()
        #else:
            #raise Exception('YAY!')
        gff3_lines.extend(align_feat_lines)
        
    # Add sort code here if needed
    #  ---- sort code here ----
    
    gff3out = open(gff3Path,'w')
    for line in gff3_lines:
        gff3out.write('%s\n' % ('\t'.join([str(x) for x in line])))
        
    gff3out.close()
    