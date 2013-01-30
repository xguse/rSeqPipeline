"""
Uses htseq-count to convert multiple sam/bam alignments to count-type data
for specific types of GTF/GFF features given an annotation.  A table file 
is created following the template below:

             aln_1    aln_2    aln_N
feature_1    int      int      int
feature_2    int      int      int
feature_3    int      int      int

"""
# TODO: => Solve the bam piped to sam problem...

import sys

import argparse

from collections import defaultdict
from collections import namedtuple

import HTSeq
#from HTSeq.scripts import count

from rSeq.utils.externals import runExternalApp
#from rSeq.scripts import my_htseq_count as count

def build_arg_str(args,sam_file):
    '''Builds argument string for external call to htseq-count based on values in args object.
    '''
    test_args = (('mode',"-m %s "),
                 ('stranded',"-s %s "),
                 ('minaqual',"-a %s"),
                 ('featuretype',"-t %s "),
                 ('idattr',"-i %s "),
                 ('samout',"-o %s "),
                 ('quiet',"-q "))
    arg_str = ''
    args_to_include = []
    
    for arg,flag in test_args:
        arg = args.__getattribute__(arg)
        if arg:
            args_to_include.append(arg)
            arg_str += flag
    # add position args to the end    
    arg_str += '%s %s'
    args_to_include.extend((sam_file,args.gff))
    
    # form final string
    arg_str = arg_str % tuple(args_to_include)
    return arg_str

def main():
    # TODO: Manage Commandline options
    argParser = argparse.ArgumentParser( 
        description=
        """This script take multiple output tables from htseq-count and combines
them into a new table with the original data as columns.""")
    

    argParser.add_argument('htseq_tables', type=str, nargs='+',
                           help='at least one htseq-count file')

    argParser.add_argument('-n',dest='table_names', type=str, nargs='+',
                           help='short names to use as table headings to represent the provided ' + 
                           'table files, otherwise file names are used, IMPORTANT: order of ' +
                           'names must match order of table files in "htseq_tables" argument!!',
                           default=False)   

    argParser.add_argument( "-o", type=str, dest="tableout",
                            default = False, help = "write out combined table to file instead of standard out")    



    args = argParser.parse_args()
    if args.table_names is not False:
        if len(args.htseq_tables) != len(args.table_names):
            raise Exception("The number of table files provided do not match the number of names provided.")
    else:
        args.table_names = args.htseq_tables

    # Map short names to SAM files
    name_map = {}
    for i,v in enumerate(args.htseq_tables):
        name_map[v] = args.table_names[i]

    # Init :out_table: (DefaultDict of Dicts.)
    out_table = defaultdict(defaultdict)

    # Iterate through table files, collect feature counts for each file and update :out_table:
    for table_file_path in args.htseq_tables:
        table_file = open(table_file_path,'rU')

        for line in table_file:
            feature,read_count = line.strip('\n').split('\t')
            out_table[feature][name_map[table_file_path]] = int(read_count)
            
        table_file.close()

    # Pull out the special features before sorting
    special_features = ('no_feature',
                        'ambiguous',
                        'too_low_aQual',
                        'not_aligned',
                        'alignment_not_unique')
    special_dict = {}
    for sf in special_features:
        special_dict[sf] = out_table.pop(sf)
        
    # Sort :out_table: keys
    sorted_features = sorted(out_table.keys())
    
    # Write out the cummulative counts
    if args.tableout is not False:
        writer = open(args.tableout,'w')
    else:
        writer = sys.stdout
    
    writer.write('feature_id\t%s\n' % ('\t'.join(args.table_names)))
    
    for feature_id in sorted_features:
        line = [feature_id]
        for table_name in args.table_names:
            line.append(str(out_table[feature_id][table_name]))
        writer.write('%s\n' % ('\t'.join(line)))
        
    for feature_id in special_features:
        line = [feature_id]
        for table_name in args.table_names:
            line.append(str(special_dict[feature_id][table_name]))
        writer.write('%s\n' % ('\t'.join(line)))  
        
    # Clean up
    if args.tableout is not False:
        writer.close()   

if __name__ == "__main__":
    main()