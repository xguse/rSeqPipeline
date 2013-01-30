"""
tool for mariangela to sort out the number of cufflinks novel
genes who have at least one exon far from any annotated exon.
"""

import argparse


def main():

    desc = """tool for mariangela to sort out the number of cufflinks novel
genes who have at least one exon far from any annotated exon."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    
    parser.add_argument('gtf', type=str,
                        help="""Path to gtf file. \n(default: %(default)s)""")
    
    args = parser.parse_args()
    
    xloc_set = set()
    
    gtf = open(args.gtf,'rU')
    
    for line in gtf:
        line = line.strip('\n').split('\t')
        comments = line[8].split('"')
        xloc = comments[1]
        xloc_set.add(xloc)
    
    print 'Number of unique xloc symbols: %s' % (len(xloc_set))
    print 'Unique xloc symbols:\n%s' % ('\n'.join(sorted(list(xloc_set))))


if __name__ == "__main__":
    main()