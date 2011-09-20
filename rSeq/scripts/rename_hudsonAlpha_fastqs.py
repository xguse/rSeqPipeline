desc = """Description: Hudson Alpha seq files come with names like this: 
(D04FKACXX_s1_1_illumina12index_1_SL6975.fastq.gz). 

We have project names for each sample like this: 
(1325-AAJ-0011). 

This script renames seq files in current directory to start with the project sample names: 
(1325-AAJ-XXXX.D04FKACXX_sX_X_illuminaXXindex_X_SLXXXX.fastq.gz). 

Using a key file formated like: 
(Index   Sequence        Library Name)
"""

import sys,os
import argparse

from rSeq.utils.files import tableFile2namedTuple

def main():
    """
    """
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('key', type=str,
                        help="""Path to file with rename data as in Description.""")
    parser.add_argument('dir', type=str,
                        help="""Path to directory with files to rename.""")
    
    
    # print the called command:
    print " ".join(sys.argv)
    
    args = parser.parse_args()
    
    keyData = tableFile2namedTuple(args.key)
    
    wrkDirFiles = os.listdir(args.dir)
    
    for f in wrkDirFiles:
        success = False
        for k in keyData:
            if k.Index and k.Library in f:
                newName = '%s.%s' % (k.Name,f)
                os.rename(f,newName)
                success = True
                print "Renamed %s as %s." % (f,newName)
                break
            else:
                pass
        if not success:
            print "%s was not in the key file." % (f)
            
if __name__ == "__main__":
    main()