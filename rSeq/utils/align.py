import sys
import subprocess
import pdb

import pysam

from rSeq.utils.errors import *
from rSeq.utils.externals import runExternalApp,mkdirp
from rSeq.utils.misc import whoami

# import pp if avail
try:
    import pp
except ImportError:
    pass


#def samCH2supCntg(samFile,outFile,cnvsnDict):
    #from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert# FIX THIS!
    #"""Convert CHxxx to superContig1.xxx"""
    
    #outFile = open(outFile, 'w')
    
    #for line in open(samFile, 'rU'):
        #if line.startswith('@SQ'):
            #line = line.split('\t')
            #line[1] = line[1].split(':')[0],':',cnvsnDict[line[1].split(':')[1]]

def cigarStr2AlignCoords(cigarString,minCoord,maxCoord,alignStrand,cigKind='ensembl',intron='I'):
    """Returns List of coordinate tuples for each block in form of:
    [
    (blkMinCoord_1,blkMaxCoord_1),
    (blkMinCoord_2,blkMaxCoord_2),
    ...,
    (blkMinCoord_n,blkMaxCoord_n)
    ]
    
    NOTE: Assumes 0-referenced "in-between" coordinate system as used in BED/exonserate."""
    
    cigInfo = parseCigarString(cigarString,kind=cigKind)
    
    # Reverse cigData if alignStrand is neg
    if alignStrand == ('-' or '-1'):
        cigInfo.reverse()
    elif alignStrand == ('+' or '1'):
        pass
    else:
        raise InvalidOptionError("cigarStr2AlignCoords: valid alignStrand values are %s. You gave: %s" \
                                 % (['-','+','-1','1'],alignStrand))
    if intron == 'I':
        ignore = 'D'
    elif intron == 'D':
        ignore = 'I'
    else:
        raise InvalidOptionError("cigarStr2AlignCoords: valid intron values are %s. You gave: %s" \
                                 % (['I','D'],intron)) 
    
    # 1: remove all 'ignore' info since these seem to relate to insertions in the aligned seq
    # and do not affect coords on the chrm seq.
    cigNoIgnores = []
    for elem in cigInfo:
        if elem[0] == ignore:
            continue
        else:
            cigNoIgnores.append(elem)
            
    # 2: combine any 'M' info that are adjacent to each other after cleaning out the 'ignore' info.
    cigCmbnMs = []
    while 1:
        if not len(cigNoIgnores) > 0:
            break
        Ms = []
        while 1:
            try:
                if cigNoIgnores[0][0] == 'M':
                    Ms.append(cigNoIgnores.pop(0))
                else:
                    break
            except IndexError:
                break
        cigCmbnMs.append(('M',sum([x[1] for x in Ms])))
        if len(cigNoIgnores) > 0:
            cigCmbnMs.append(cigNoIgnores.pop(0))
        else:
            break
    
    # 3: walk the cig, keeping track of where on Chrm you are and logging boundaries.
    blkCoordsRelative = []
    i = 0
    while i < len(cigCmbnMs):
        if cigCmbnMs[i][0] == 'M':
            end   = sum([x[1] for x in cigCmbnMs[:i+1]]) 
            start = end - cigCmbnMs[i][1]
            blkCoordsRelative.append((start,end))
            i += 1
        else:
            i += 1
    blkCoordsAbsolute = []
    for coord in blkCoordsRelative:
        blkCoordsAbsolute.append((coord[0]+minCoord,coord[1]+minCoord))
    
    # Confirm that calculated max coord == maxCoord
    calcMax = max([x[1] for x in blkCoordsAbsolute])
    if calcMax != maxCoord:
        raise UnexpectedValueError("cigarStr2AlignCoords: calculated max coord != maxCoord.  My code may not be correct, but check that you used the correct 'intron' value first.")
    else:
        return blkCoordsAbsolute
    
def parseCigarString(cigarString,kind='ensembl'):
    """Returns cigarInfo from cigarString in form:
    [
    [operationLetter_1,run_1],
    ...,
    [operationLetter_n,run_n]
    ]
    
    Validates operation letters before returning."""
    
    validKinds  = ['ensembl','exonerate']
    validCigOps = ['M','D','I']
    if kind not in validKinds:
        raise InvalidOptionError('parseCigarString: valid kinds are: %s.  You gave: %s' % (validKinds,kind))
    if not type(cigarString) == type(''):
        raise InvalidOptionError('parseCigarString: type(cigarString) != type(""): %s' % (kind))
    
    cigData = None
    if kind == validKinds[0]:
        cigData = parseEnsemblCigar(cigarString)
    else:
        cigData = parseExonerateCigar(cigarString)
    
    # Do some validation of cigar operation letters:
    cigOps = list(set([x[0] for x in cigData]))
    for op in cigOps:
        if not op in validCigOps:
            raise InvalidOptionError("parseCigarString: valid cigar operation letters are %s.  malformed cigarString (%s)." \
                                     % (validCigOps,cigarString))
    return cigData
    

def parseExonerateCigar(cigarString):
    """Given exonerate style cigar string, return list of lists representing
    [[letter,number],[letter,number],...]."""
    cigarList = cigarString.split()
    
    # Format cigar info for easy use:
    cigar = []
    if not len(cigarList)%2==0:
        raise InvalidOptionError('parseExonerateCigar: cigarList var is not a list with an even number of indexes.')
    while cigarList:
        cigar.append((cigarList.pop(0),int(cigarList.pop(0))))
    return cigar

def parseEnsemblCigar(cigarString):
    """Given ensembl style cigar line, return list of lists representing
    [[letter,number],[letter,number],...]."""
    def getNumEnd(cigString,index):
        try:
            int(cigString[index])
            index+=1
            end = getNumEnd(cigString,index)
        except ValueError:
            return index
        return end
    
    cigList = []
    i=0
    while i < len(cigarString):
        num = None
        let = None
        
        try:
            int(cigarString[i])
            try:
                end = getNumEnd(cigarString,i)
                num = int(cigarString[i:end])
                i = end
            except ValueError:
                raise Exception('parseEnsemblCigar.getNumEnd: failed check code.')
        except ValueError:
            let = cigarString[i]
            i+=1
            
        if num:
            let = cigarString[i]
            i+=1
            cigList.append((let,num))
        elif let:
            cigList.append((let,1))
        else:
            raise InvalidOptionError('parseEnsemblCigar: passed through a cycle without assigning let OR num. Check for valid cigarString: %s' % (cigarString))
    return cigList
            

def exonerateCigar2BEDline(cigarLine,rgb="0,0,0"):
    """Takes a cigar line from an exonerate output file and returns
    a BED formated line representing the alignment blocks of the query
    as aligned to the target. """
    
    line  = cigarLine.strip('\n').split() # its a space delimited line
    
    # exonerate data
    query   = line[1]   # query id
    qStart  = line[2]   # ! will be lower number than qEnd if qStrand is '-'
    qEnd    = line[3]   # ! will be higher number than qStart if qStrand is '-'
    qStrand = line[4]   #
    target  = line[5]   # target id
    tStart  = line[6]   # ! will be lower number than tEnd if tStrand is '-'
    tEnd    = line[7]   # ! will be higher number than tStart if tStrand is '-'
    tStrand = line[8]   #
    score   = line[9]   # "sum of transistion scores used in the dynamic programming."
    cigInfo = ' '.join(line[10:])  # actual CIGAR encoding converted back to string for parser.
    
    # Format cigar info for easy use:
    cigar = parseCigarString(cigInfo,type='exonerate')
    if "I" in [x[0] for x in cigar]:
        sys.stderr.write('warn: just encountered a cigar containing "I".  Skipping.\n')
        return '' # for now just do nothing

        #raise UnexpectedValueError('cigar info includeds "I" and I have not been taught how to deal with this!')
    
    # If match is to '-' strand of target, reverse the cigar info
    if tStrand == '-':
        cigar.reverse()
        
    # Bed Terms
    chrom       = target
    chromStart  = min(int(tStart),int(tEnd))
    chromEnd    = max(int(tStart),int(tEnd))
    name        = query
    score       = 0
    strand      = tStrand
    thickStart  = chromStart
    thickEnd    = chromEnd
    itemRgb     = rgb
    blockCount  = len([x[0] for x in cigar if x[0]=="M"])
    blockSizes  = ','.join([x[1] for x in cigar if x[0]=="M"])
    blockStarts = ['0']  # for now.  See below.
    
    # Calculate blockStarts
    i = 0
    while 1:
        if i < len(cigar)-1:
            blockStarts.append(str(int(cigar[i][1])+int(cigar[i+1][1])))
            i += 2
        else:
            break
    blockStarts = ','.join(blockStarts)
    
    # Formulate and return BedLine
    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
           (chrom,chromStart,chromEnd,name,score,strand,
            thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)
    

############################
########## BOWTIE ##########
############################

def pipe_bowtie2bam(btIndex,outPath,fastq1,fastq2=None,options=None):
    """Run bowtie and pipe output to samtools to convert to
    BAM type.
    
    btIndex: valid bowtie index residing in $BOWTIE_INDEXES
    outPath: path to write the output
    fastq1: path(s) to upstream mates or single reads (if fastq2=None) ["p1,p2...pN"]
    fastq2: path(s) to downstream mates ["p1,p2...pN"]
    options: valid cmd line options string for bowtie 
    
    Exp of options (mimic ELAND and write as SAM):
    options = "--solexa1.3-quals -v 2 -m 1 -S"
    
    ----------
    bowtie help text:
    
    Usage: 
    bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
  
    <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
    <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
    <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-" for stdin.
    <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-" for stdin.
    <hit>   File to write hits to (default: stdout)
    Input:
      -q                 query input files are FASTQ .fq/.fastq (default)
      -f                 query input files are (multi-)FASTA .fa/.mfa
      -r                 query input files are raw one-sequence-per-line
      -c                 query sequences given on cmd line (as <mates>, <singles>)
      -C                 reads and index are in colorspace
      -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
      --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
      -s/--skip <int>    skip the first <int> reads/pairs in the input
      -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
      -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
      -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
      --phred33-quals    input quals are Phred+33 (default)
      --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
      --solexa-quals     input quals are from GA Pipeline ver. < 1.3
      --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
      --integer-quals    qualities are given as space-separated integers (not ASCII)
    Alignment:
      -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
        or
      -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
      -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
      -l/--seedlen <int> seed length for -n (default: 28)
      --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
      -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
      -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
      --nofw/--norc      do not align to forward/reverse-complement reference strand
      --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
      --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
      -y/--tryhard       try hard to find valid alignments, at the expense of speed
      --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
    Reporting:
      -k <int>           report up to <int> good alignments per read (default: 1)
      -a/--all           report all alignments per read (much slower than low -k)
      -m <int>           suppress all alignments if > <int> exist (def: no limit)
      -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
      --best             hits guaranteed best stratum; ties broken by quality
      --strata           hits in sub-optimal strata aren't reported (requires --best)
    Output:
      -t/--time          print wall-clock time taken by search phases
      -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
      --quiet            print nothing but the alignments
      --refout           write alignments to files refXXXXX.map, 1 map per reference
      --refidx           refer to ref. seqs by 0-based index rather than name
      --al <fname>       write aligned reads/pairs to file(s) <fname>
      --un <fname>       write unaligned reads/pairs to file(s) <fname>
      --max <fname>      write reads/pairs over -m limit to file(s) <fname>
      --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
      --fullref          write entire ref name (default: only up to 1st space)
    Colorspace:
      --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
         or
      --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
      --col-cseq         print aligned colorspace seqs as colors, not decoded bases
      --col-cqual        print original colorspace quals, not decoded quals
      --col-keepends     keep nucleotides at extreme ends of decoded alignment
    SAM:
      -S/--sam           write hits in SAM format
      --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
      --sam-nohead       supppress header lines (starting with @) for SAM output
      --sam-nosq         supppress @SQ header lines for SAM output
      --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
    Performance:
      -o/--offrate <int> override offrate of index; must be >= index's offrate
      -p/--threads <int> number of alignment threads to launch (default: 1)
      --mm               use memory-mapped I/O for index; many 'bowtie's can share
      --shmem            use shared mem for index; many 'bowtie's can share
    Other:
      --seed <int>       seed for random number generator
      --verbose          verbose output (for debugging)
      --version          print version information and quit
      -h/--help          print this usage message

    """
    Popen = subprocess.Popen
    PIPE = subprocess.PIPE
    
    # example:
    # bowtie --solexa1.3-quals -v 2 -m 1 -S $BINDX -p 2 --mm "${RBA},${RBB}" | samtools view -bS - > $RB_BAM
    
    # Prepare the command line args
    options = options.split()
    #pdb.set_trace()
    if fastq2 == None:
        btArgs = ['bowtie'] + options + [btIndex] + [fastq1]
        samArgs = ['samtools','view','-bS','-']
        print "running: \nbowtie with cmd: '%s'\nsamtools with cmd: '%s'\noutfile:%s" \
              % (' '.join(btArgs),' '.join(samArgs),outPath)
        
        bowtieStep  = Popen(btArgs,
                            stdout=PIPE)
        
        samViewStep = Popen(samArgs,
                            stdin=bowtieStep.stdout,
                            stdout=open(outPath,'w'))

        bowtieStep.stdout.close()
        samViewStep.stdout.close()
    else:
        btArgs = ['bowtie'] + options + [btIndex] + ['-1'] + [fastq1] + ['-2'] + [fastq2]
        samArgs = ['samtools','view','-bS','-']
        print "running: \nbowtie with cmd: '%s'\nsamtools with cmd: '%s'\noutfile:%s" \
              % (' '.join(btArgs),' '.join(samArgs),outPath)
        
        bowtieStep  = Popen(args,
                            stdout=PIPE)
        
        samViewStep = Popen(samArgs,
                            stdin=bowtieStep.stdout,
                            stdout=open(outPath,'w'))

        bowtieStep.stdout.close()
        samViewStep.stdout.close()
        
    

def bowtie_index(reference_in,ebwt_outfile_base,runDir,options=None):
    """Create bowtie indexes from new fasta set.
    options : quoted string representing valid cmd line bowtie-build options
    runDir  : path to dir to place stdErr/stdOut logs - all steps of pipeline scripts should share same runDir
    
    ----------
    bowtie-build help text:
    
    Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>
        reference_in            comma-separated list of files with ref sequences
        ebwt_outfile_base       write Ebwt data to files with this dir/basename
    Options:
        -f                      reference files are Fasta (default)
        -c                      reference sequences given on cmd line (as <seq_in>)
        -C/--color              build a colorspace index
        -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
        -p/--packed             use packed strings internally; slower, uses less mem
        -B                      build both letter- and colorspace indexes
        --bmax <int>            max bucket sz for blockwise suffix-array builder
        --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
        --dcv <int>             diff-cover period for blockwise (default: 1024)
        --nodc                  disable diff-cover (algorithm becomes quadratic)
        -r/--noref              don't build .3/.4.ebwt (packed reference) portion
        -3/--justref            just build .3/.4.ebwt (packed reference) portion
        -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)
        -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
        --ntoa                  convert Ns in reference to As
        --seed <int>            seed for random number generator
        -q/--quiet              verbose output (for debugging)
        -h/--help               print detailed description of tool and its options
        --usage                 print this usage message
        --version               print version information and quit

    """
    
    # make runDir if it does not yet exist
    print "creating: %s" % (runDir)
    mkdirp(runDir)
    
    # Construct cmdArgs
    if options:
        cmdArgs = "%s %s %s" % (options,reference_in,ebwt_outfile_base)
    # Run bowtie-build and capture output
    btBuildResults = runExternalApp('bowtie-build',cmdArgs)
    
    # Report bowtie-build stdout and stderr
    if btBuildResults[0]:
        print('[%s] %s' % (whoami(),btBuildResults[0]))
    if btBuildResults[1]:
        sys.stderr('[%s] %s' % (whoami(),btBuildResults[1]))
        
    return btBuildResults
    
def bowtie_align(ebwt,readsString,hit,runDir,options=None):
    """Run alignment of fastQ to bowtie index.
    options     : quoted string representing valid cmd line bowtie-build options
    runDir      : path to dir to place stdErr/stdOut logs - all steps of pipeline scripts should share same runDir
    readsString : appropriate quoted string representing which fastq files to use (see bowtie -h).
    
    ----------
    bowtie help text:
    
    Usage: 
    bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
  
    <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
    <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
    <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-" for stdin.
    <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-" for stdin.
    <hit>   File to write hits to (default: stdout)
    Input:
      -q                 query input files are FASTQ .fq/.fastq (default)
      -f                 query input files are (multi-)FASTA .fa/.mfa
      -r                 query input files are raw one-sequence-per-line
      -c                 query sequences given on cmd line (as <mates>, <singles>)
      -C                 reads and index are in colorspace
      -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
      --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
      -s/--skip <int>    skip the first <int> reads/pairs in the input
      -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
      -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
      -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
      --phred33-quals    input quals are Phred+33 (default)
      --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
      --solexa-quals     input quals are from GA Pipeline ver. < 1.3
      --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
      --integer-quals    qualities are given as space-separated integers (not ASCII)
    Alignment:
      -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
        or
      -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
      -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
      -l/--seedlen <int> seed length for -n (default: 28)
      --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
      -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
      -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
      --nofw/--norc      do not align to forward/reverse-complement reference strand
      --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
      --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
      -y/--tryhard       try hard to find valid alignments, at the expense of speed
      --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
    Reporting:
      -k <int>           report up to <int> good alignments per read (default: 1)
      -a/--all           report all alignments per read (much slower than low -k)
      -m <int>           suppress all alignments if > <int> exist (def: no limit)
      -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
      --best             hits guaranteed best stratum; ties broken by quality
      --strata           hits in sub-optimal strata aren't reported (requires --best)
    Output:
      -t/--time          print wall-clock time taken by search phases
      -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
      --quiet            print nothing but the alignments
      --refout           write alignments to files refXXXXX.map, 1 map per reference
      --refidx           refer to ref. seqs by 0-based index rather than name
      --al <fname>       write aligned reads/pairs to file(s) <fname>
      --un <fname>       write unaligned reads/pairs to file(s) <fname>
      --max <fname>      write reads/pairs over -m limit to file(s) <fname>
      --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
      --fullref          write entire ref name (default: only up to 1st space)
    Colorspace:
      --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
         or
      --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
      --col-cseq         print aligned colorspace seqs as colors, not decoded bases
      --col-cqual        print original colorspace quals, not decoded quals
      --col-keepends     keep nucleotides at extreme ends of decoded alignment
    SAM:
      -S/--sam           write hits in SAM format
      --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
      --sam-nohead       supppress header lines (starting with @) for SAM output
      --sam-nosq         supppress @SQ header lines for SAM output
      --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
    Performance:
      -o/--offrate <int> override offrate of index; must be >= index's offrate
      -p/--threads <int> number of alignment threads to launch (default: 1)
      --mm               use memory-mapped I/O for index; many 'bowtie's can share
      --shmem            use shared mem for index; many 'bowtie's can share
    Other:
      --seed <int>       seed for random number generator
      --verbose          verbose output (for debugging)
      --version          print version information and quit
      -h/--help          print this usage message

    """
    # make runDir if it does not yet exist
    mkdirp(runDir)
    
    # Construct cmdArgs
    if options:
        cmdArgs = "%s %s  %s %s" % (options,ebwt,readsString,hit)
    else:
        cmdArgs = "%s %s %s" % (ebwt,readsString,hit)
        
    # Run and capture output
    print "Setting up bowtie call with the following cmd:\n\t\tbowtie %s" % (cmdArgs)
    btResults = runExternalApp('bowtie',cmdArgs)
    
    # Report stdout and stderr
    if btResults[0]:
        for line in btResults[0].split('\n'):
            print('[%s] %s' % (whoami(),line))
        
    if btResults[1]:
        for line in btResults[1].split('\n'):
            sys.stderr.write('[%s] %s' % (whoami(),line))
        
    return btResults

##############################
########## SAMTOOLS ##########
##############################

def run_samtools(tool,argsList):
    """Core wrapper for calls to samtools and associated
    stderr/stdout reporting.  See specific func defs for each
    tool's samtool's help text.
    
    tool     : sort,view,index,etc...
    argsList : list obj with appropriatly structured args for samtools "tool"
    """
    # Construct cmdArgs
    cmdArgs = "%s %s" % (tool," ".join(argsList))
        
    # Run and capture output
    print "Setting up samtools call with the following cmd:\n\tsamtools %s" % (cmdArgs)
    results = runExternalApp('samtools',cmdArgs)
    
    # Report stdout and stderr
    if results[0]:
        for line in results[0].split('\n'):
            print('[%s] %s' % (whoami(),line))
        
    if results[1]:
        for line in results[1].split('\n'):
            sys.stderr.write('[%s] %s' % (whoami(),line))
        
    return results

def samtools_view(argsList):
    """Wrapper for samtools "view".
    
    argsList  :  list obj with appropriatly structured args for samtools "view".
    
    ----------
    samtools view help text:
    
    Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]

    Options: -b       output BAM
             -h       print header for the SAM output
             -H       print header only (no alignments)
             -S       input is SAM
             -u       uncompressed BAM output (force -b)
             -x       output FLAG in HEX (samtools-C specific)
             -X       output FLAG in string (samtools-C specific)
             -c       print only the count of matching records
             -t FILE  list of reference names and lengths (force -S) [null]
             -T FILE  reference sequence file (force -S) [null]
             -o FILE  output file name [stdout]
             -R FILE  list of read groups to be outputted [null]
             -f INT   required flag, 0 for unset [0]
             -F INT   filtering flag, 0 for unset [0]
             -q INT   minimum mapping quality [0]
             -l STR   only output reads in library STR [null]
             -r STR   only output reads in read group STR [null]
             -?       longer help

    
    """
        
    return run_samtools('view',argsList)

def samtools_sort(argsList):
    """Wrapper for samtools "sort".
    
    argsList  :  list obj with appropriatly structured args for samtools "sort".
    
    ----------
    samtools sort help text:
    
    Usage: samtools sort [-on] [-m <maxMem>] <in.bam> <out.prefix>
    """
    
    return run_samtools('sort',argsList)

def samtools_index(argsList):
    """Wrapper for samtools "index".
    
    argsList  :  list obj with appropriatly structured args for samtools "index".
    
    ----------
    samtools sort help text:
    
    Usage: samtools index <in.bam> [out.index]
    """
    
    return run_samtools('index',argsList)

def pysamtools_call(funcName,argsList):
    """
    ALERT: pysam's encapsulation of AT LEAST 'samtools view'
        hijack's stderr/stdout which among other problems (unfaithful interpretation of commandline 
        arguments like view's "-o" and "-S") make the following code unreliable for now.
    
    
    Wrapper for pysam."funcName"() that adds useful reporting information.
    
    argsList  :  list obj with appropriatly structured args for pysam."funcName"()
    """
    if 1==1:
        raise NotImplementedError("""ERROR: pysam's encapsulation of AT LEAST 'samtools view'
        hijack's stderr/stdout which among other problems (unfaithful interpretation of commandline 
        arguments like view's "-o" and "-S") make the following code unreliable for now.""")
    # Make sure that the funcName exists in samtools
    if funcName not in pysam.SAMTOOLS_DISPATCH:
        raise InvalidOptionError(optVal=funcName,optName="funcName",validVals=pysam.SAMTOOLS_DISPATCH.keys())
        
    print "Setting up samtools call with the following cmd:\n\tsamtools %s %s" % (funcName,' '.join(argsList)) 
    # Construct correct pysam tool agent type and run the command
    dispatcher = pysam.SAMTOOLS_DISPATCH[funcName]
    pysamAgent = pysam.SamtoolsDispatcher(dispatcher[0],dispatcher[1])
    
    results = pysamAgent(' '.join(['"%s"' % (x) for x in argsList]))

    # Communicate StdOut
    for line in results[2].split('\n'):
        print('[%s:%s] %s' % (whoami(),funcName,line))
    # Communicate StdErr
    for line in results[1].split('\n'):
        sys.stderr.write('[%s:%s] %s' % (whoami(),funcName,line))
    
    return results


############################
########## TOPHAT ##########
############################

def tophat_align(bowtie_index,readsA,readsB=None,qualsA=None,qualsB=None,options=None,runDir=None):
    # TODO: eliminate "Reconstituting reference FASTA file from Bowtie index \ [FAILED] \ Error: bowtie-inspect returned an error."

    """
    Wrapper for calling tophat and dealing with output.
    
    **ARGS**
    _name_         _type_      _desc_
    bowtie_index : String      bowtie index base-name in $BOWTIE_INDEXES
    readsA       : List        List of FilePaths to fastQ files
    readsB       : List/None   
    qualsA       : List/None
    qualsB       : List/None
    options      : String      quoted comma-sep str of CLI tophat options
    ------------------
    
    **NOTES**
    If (readsA AND readsB) OR (qualsA AND qualsB), the order of filePaths
    corresponding to PE-mates must match:
        readsA,readsB = [readsA_1,readsA_2], [readsB_1,readsB_2]
        qualsA,qualsB = [qualsA_1,qualsA_2], [qualsB_1,qualsB_2]
    ------------------
    
    
    **TOPHAT HELP TEXT FOLLOWS**
tophat: 
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat [options] <bowtie_index> <reads1[,reads2,...,readsN]> [reads1[,reads2,...,readsN]] [quals1,[quals2,...,qualsN]] [quals1[,quals2,...,qualsN]]
    
Options:
    -v/--version
    -o/--output-dir                <string>    [ default: ./tophat_out ]
    -a/--min-anchor                <int>       [ default: 8            ]
    -m/--splice-mismatches         <0-2>       [ default: 0            ]
    -i/--min-intron-length         <int>       [ default: 50           ]
    -I/--max-intron-length         <int>       [ default: 500000       ]
    -g/--max-multihits             <int>       [ default: 40           ]
    -F/--min-isoform-fraction      <float>     [ default: 0.15         ]
    --max-insertion-length         <int>       [ default: 3            ]
    --max-deletion-length          <int>       [ default: 3            ]
    --solexa-quals                          
    --solexa1.3-quals                          (same as phred64-quals)
    --phred64-quals                            (same as solexa1.3-quals)
    -Q/--quals
    --integer-quals
    -C/--color                                 (Solid - color space)
    --color-out
    --library-type                             (--fr-unstranded, --fr-firststrand, --fr-secondstrand, --ff-unstranded, --ff-firststrand, --ff-secondstrand)
    -p/--num-threads               <int>       [ default: 1            ]
    -G/--GTF                       <filename>
    -j/--raw-juncs                 <filename>
    --insertions		   <filename>
    --deletions			   <filename>
    -r/--mate-inner-dist           <int>       
    --mate-std-dev                 <int>       [ default: 20           ]
    --no-novel-juncs
    --allow-indels
    --no-novel-indels                           
    --no-gtf-juncs                             
    --no-coverage-search
    --coverage-search                                              
    --no-closure-search
    --closure-search
    --fill-gaps        
    --microexon-search
    --butterfly-search
    --no-butterfly-search
    --keep-tmp
    --tmp-dir                      <dirname>
    
Advanced Options:

    --segment-mismatches           <int>       [ default: 2            ]
    --segment-length               <int>       [ default: 25           ]
    --min-closure-exon             <int>       [ default: 100          ]
    --min-closure-intron           <int>       [ default: 50           ]
    --max-closure-intron           <int>       [ default: 5000         ]
    --min-coverage-intron          <int>       [ default: 50           ]
    --max-coverage-intron          <int>       [ default: 20000        ]
    --min-segment-intron           <int>       [ default: 50           ]
    --max-segment-intron           <int>       [ default: 500000       ] 

SAM Header Options (for embedding sequencing run metadata in output):
    --rg-id                        <string>    (read group ID)
    --rg-sample                    <string>    (sample ID)
    --rg-library                   <string>    (library ID)
    --rg-description               <string>    (descriptive string, no tabs allowed)
    --rg-platform-unit             <string>    (e.g Illumina lane ID)
    --rg-center                    <string>    (sequencing center name)
    --rg-date                      <string>    (ISO 8601 date of the sequencing run)
    --rg-platform                  <string>    (Sequencing platform descriptor)  

    for detailed help see http://tophat.cbcb.umd.edu/manual.html

    """
    
    # make runDir if it does not yet exist
    #   => tophat takes care of this for us
    ##mkdirp(runDir)
    
    # Construct cmdArgs
    if not type(readsA) == type([]):
        raise TypeError('readsA type should be "[]".')
    readsA = ','.join(readsA)
    
    # format readsB
    if readsB:
        if not type(readsB) == type([]):
            raise TypeError('readsB type should be "[]".')
        readsB = ','.join(readsB)
    else:
        readsB = ''
    
    # format qualsA    
    if qualsA:
        if not type(qualsA) == type([]):
            raise TypeError('qualsA type should be "[]".')
        qualsA = ','.join(qualsA)
    else:
        qualsA = ''
    
    # format qualsB
    if qualsB:
        if not type(qualsB) == type([]):
            raise TypeError('qualsB type should be "[]".')
        qualsB = ','.join(qualsB)
    else:
        qualsB = ''
    
    # format runDir
    if runDir:
        runDir = ' -o %s' % (runDir)
    else:
        runDir = ''
    
    # format cmdArgs
    if options == None:
        options = ''
        
    # format bowtie_index
    #   This seems needed or tophat cant reconstitute the fasta files for some reason
    #   and returns the following error:
    #    *  Reconstituting reference FASTA file from Bowtie index
    #    *  [FAILED]
    #    *  Error: bowtie-inspect returned an error.
    if not bowtie_index.startswith('/'):
        bowtie_index = "$BOWTIE_INDEXES/%s" 
        
    cmdArgs = "%s %s %s %s %s %s" % (options+runDir,bowtie_index,readsA,readsB,qualsA,qualsB)

        
    # Run and capture output
    print "Setting up tophat call with the following cmd:\n\t\ttophat %s" % (cmdArgs)
    thResults = runExternalApp('tophat',cmdArgs)
    
    # Report stdout and stderr
    if thResults[0]:
        for line in thResults[0].split('\n'):
            print('[%s] %s' % (whoami(),line))
        
    if thResults[1]:
        for line in thResults[1].split('\n'):
            sys.stderr.write('[%s] %s' % (whoami(),line))
        
    return thResults

