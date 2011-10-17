import sys,os
import argparse
import multiprocessing as mp


from rSeq.utils.errors import *
from rSeq.utils.files import filter_PEfastQs
from rSeq.utils.externals import mkdirp

def main():
    """Inputs:
    -- Txt table file containing some tab-delim file path inputs for files.filter_PEfastQs()
           PE_FastqPathFwd<tab>PE_FastqPathRev<tab>OutputComboBaseName<newLine>
    -- String representing a lambda func to act as filter for fastqRecs
    Outputs:
    -- Writes filtered data to paths specified in the input txt table file"""
    
    desc  = """This script filters paired fastq files based on a provided lambda filter logic.  Input and output paths are determined by the input table file."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('input_table', type=str,
                        help="""Path to input table file.""")
    parser.add_argument('filter_func', type=str,
                        help="""lambda filter function (quoted string).""")
    parser.add_argument('out_dir', type=str,
                        help="""Path to out directory.""")





    
    
    args = parser.parse_args()
    
    # print the called command:
    sys.stderr.write("%s\n" % (" ".join(sys.argv)))
    
    # open and parse the input table file
    inputs = [x.strip('\n').split('\t') for x in open(args.input_table,'rU')]
    for i in inputs:
        if len(i) != 3:
            raise InvalidFileFormatError("At least one line in %s does not have exactly three columns:\n%s\n" % (args.input_table,i))
    
    # create out_dir if needed
    mkdirp(args.out_dir)
    
    # set up and unleash the subprocesses
    filtFunc = eval(args.filter_func,{"__builtins__":None})
    jobs = []
    for i in inputs:
        arguments = [filtFunc, i[0], i[1]]
        # build the output file names (matchedPassPath1,matchedPassPath2,singlePassPath,nonPassPath) from baseNames 
        arguments.append("%s/%s.filtered.mated.fastq"   % (args.out_dir.rstrip('/'),i[0].split('/')[-1].replace('.fastq','')))   # matchedPassPath1
        arguments.append("%s/%s.filtered.mated.fastq"   % (args.out_dir.rstrip('/'),i[1].split('/')[-1].replace('.fastq','')))   # matchedPassPath2
        arguments.append("%s/%s.filtered.singled.fastq" % (args.out_dir.rstrip('/'),i[2].split('/')[-1].replace('.fastq',''))) # singlePassPath
        arguments.append("%s/%s.filtered.failed.fastq"  % (args.out_dir.rstrip('/'),i[2].split('/')[-1].replace('.fastq','')))  # nonPassPath
        
        p = mp.Process(target=filter_PEfastQs,args=tuple(arguments))
        jobs.append(p)
        p.start()
    
    return jobs
    
    
    


        
if __name__ == "__main__":
    jobs = main()
    # Require that all jobs be completed before going on
    for j in jobs:
        j.join()
    # check to ensure that none of the subprocesses returned a non-zero exit status
    for j in jobs:
        if j.exitcode != 0:
            raise Exception("Job %s returned a non-zero exit status: %s" % (j.name,j.exitcode))