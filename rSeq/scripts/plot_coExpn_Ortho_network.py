"""
1: ...
"""

import sys,os
import argparse
import logging
import math

from matplotlib import pyplot as plt
# Something is going wrong if I just try to import nx
# This fixes it. Don't ask me what is going on...
try:
    from xml.etree import cElementTree
except AttributeError:
    pass
import networkx as nx
from networkx import graphviz_layout

from rSeq.utils.errors import *
from rSeq.utils.files import tableFile2namedTuple
from rSeq.utils.networks import import_edges
from rSeq.utils.networks import graph_connected_nodes
from rSeq.utils.networks import weight_edges_with_pearsonr
from rSeq.utils.expression import mangle_expn_vectors



def main():
    """
    1: ...
    """
    
    desc = """... ask me later! I'm on a deadline! ..."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    logger = logging.getLogger(sys.argv[0].split('/')[-1])
    
    parser.add_argument('--hfile', type=str, required=True, nargs='+',
                        help="""Quoted ;-delimited list containing info about the files containing homology relationships:
[--hfile "path;header1;header2"].  At LEAST one hfile is required and all MUST have all three trailing data.""")
    
    parser.add_argument('--xprnfile', type=str, required=True, nargs='+',
                        help="""Quoted ;-delimited list containing info about the files containing expression data:
[--xprnfile "path;nameHeader;conditionHeader1;...;conditionHeaderN"].  At LEAST one xprnfile is required and all MUST have 
exactly one <path>, exactly one <nameHeader>, and at LEAST one <conditionHeader>. It is VERY important that 
you list the same number of conditions for each expnfile and that the order reflects which condition values are to be compared.""")
    
    parser.add_argument('--targets', type=str, required=True, nargs='+',
                        help="""A list of the gene/tx/protein symbols to use for pulling out all connected nodes.""")
    
    #parser.add_argument('--pfixlen', type=str, required=True, nargs='+',
                        #help="""One length of the symbol prefixes (AAEL for AAEL007639-PA) or a list of prefix le.""")
    
    parser.add_argument('--cmap', type=str, required=True, nargs='+',
                        help="""A list of species-prefix:color combinations to set the node colors:
[--cmap <AAEL:b ...>]. the number of combinations should match the number of files given to --xprnfile.""")
    
    parser.add_argument('--out', type=str, required=True,
                        help="""Path to outfile.  Its file extention chooses the file type.""")
    
    parser.add_argument('--graphml', type=str, required=False,
                        help="""Include a file path if you would like a graphML version of the final graph. (optional)""")
    
    parser.add_argument('--nonames', action='store_true',
                        help="""If used: gene/tx names will NOT be displayed.""")
    
    parser.add_argument('--noshow', action='store_true',
                        help="""If used: the graph NOT be displayed interactively.""")
    
    args = parser.parse_args()
    
    # some manual arg set-up and checking
    for i in range(len(args.hfile)):
        args.hfile[i] = args.hfile[i].split(';')
        if len(args.hfile[i]) != 3:
            raise SanityCheckError('EXACTLY 3 values must follow --hfile: you gave %s' % (args.hfile[i]))
        
    xLen = set()
    for i in range(len(args.xprnfile)):
        args.xprnfile[i] = args.xprnfile[i].split(';')
        if not len(args.xprnfile[i]) >= 3:
            raise SanityCheckError('At LEAST 3 values must follow --xprnfile: you gave %s' % (args.xprnfile[i]))
        else:
            xLen.add(len(args.xprnfile[i]))
    if not len(xLen) == 1:
        raise SanityCheckError('The same number of values must follow every --xprnfile flag.')
    
    if not len(args.xprnfile) == len(args.cmap):
        raise SanityCheckError('The length of values following --xprnfile and --cmap must be the same.')
    
    cDict = {}
    for combo in args.cmap:
        try:
            prefix,color = combo.split(':')
        except:
            raise
        cDict[prefix] = color
    
    
    
    # read in the expression vector data
    tmpDict = {}
    xDict = {}
    for xfile in args.xprnfile:
        tmpDict.update(mangle_expn_vectors(expnPath=xfile[0],txNameHeader=xfile[1],condHeaders=xfile[2:],manualHeaders=False))
        
    # convert -RX into -PX
    for k,v in tmpDict.iteritems():
        xDict[k.replace('-R','-P')] = v
        
    del(tmpDict)
    
    # lets get started: init the graph
    graph = nx.Graph()
    
    for f in args.hfile:
        import_edges(graphObj=graph,edgeTablePath=f[0],startNodeHeader=f[1],endNodeHeader=f[2])
    
        
    # remove the '' node caused by unpaired relationships
    graph.remove_node('')
    
    # for debugging
    nx.write_gpickle(graph,"/tmp/ortho1.gpickle")    
        
    # Cut out a subgraph using provided targets
    subgraph = graph_connected_nodes(graphObj=graph,nodeList=args.targets)
    
    
    # weight the edges in subgraph by the pearsonr between their expression vectors
    weight_edges_with_pearsonr(graphObj=subgraph,dataVectors=xDict,uni=False)
    
    # if the edge length is imposible to graph (inf or nan) kill the edge
    badEdges = []
    edgesMissingNodes = []
    for i,j in subgraph.edges_iter():
        try:
            if math.isnan(subgraph[i][j]['rVal']) or math.isinf(subgraph[i][j]['rVal']):
                badEdges.append((i,j))
        except KeyError:
            edgesMissingNodes.append((i,j))
                
    subgraph.remove_edges_from(badEdges)
    subgraph.remove_edges_from(edgesMissingNodes)
    
    
    # begin drawing the graph by setting the node positions
    #a = nx.to_agraph(subgraph)
    #a.layout()
    #a.draw('%s.gv.png' % args.out)
    #h = nx.from_agraph(a)
    #pos = nx.graphviz_layout(h)
    #pos= nx.spring_layout(subgraph,iterations=100)
    pos = nx.graphviz_layout(subgraph, args='-LC1000000000')
    
    # set node colors
    nodelist = subgraph.nodes()
    node_colors = []
    prefixes = cDict.keys()
    
    # get edge lables:
    eLab = {}
    for i,j in subgraph.edges_iter():
        eLab[(i,j)] = 'r = %s\np = %s' % (round(subgraph[i][j]['rVal'],3),round(subgraph[i][j]['pVal'],3))
    
    for n in nodelist:
        # crazy list comprehension python-voodoo to create a list of colors in the same order as nodelist
        node_colors.extend([cDict[x] for x in prefixes if n.startswith(x)])
    nx.draw_networkx_nodes(subgraph, pos,nodelist,node_color=node_colors,node_size=1000,node_shape='o', aplha=.7)
    sigEdges = []
    nonSigEdges = []
    for e in subgraph.edges_iter():
        if float(subgraph[e[0]][e[1]]['pVal']) <= 0.05:
            sigEdges.append(e)
        else:
            nonSigEdges.append(e)
            
    # Define color map for edge 'heats'
    g2r = {'green': ((0.0, 0.0, 0.0),
                     (0.66, 0.0, 0.0),
                     (1.0, 1.0, 1.0)),
 
          'blue': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
 
          'red':  ((0.0, 1.0, 1.0),
                   (0.33, 0.0, 0.0),
                   (1.0, 0.0, 0.0))} 
    
    
    b2g2y2o2r = {'red':   ((0.0,  0.0, 0.0),
                          #(0.9,  1.0, 1.0),
                          (1.0,  1.0, 1.0)),
       
                'green': ((0.0,  0.0, 0.0),
                          (0.4, 1.0, 1.0),
                          (0.6, 1.0, 1.0),
                          (1.0, 0.0, 0.0)),
       
                'blue':  ((0.0,  1.0,1.0),
                          #(0.1,  0.0, 0.0),
                          (1.0,  0.0, 0.0))}
    
    plt.register_cmap(name='corrMap', data=b2g2y2o2r)
    corrMap = plt.get_cmap('corrMap')
    
            
    
    nx.draw_networkx_edges(subgraph, pos, edgelist=nonSigEdges, width=2.0, edge_cmap=corrMap,
                           edge_vmin=-1,
                           edge_vmax=1,                           
                           edge_color=[subgraph[e[0]][e[1]]['weight'] for e in nonSigEdges],
                           style='dashed', alpha=.7)
    nx.draw_networkx_edges(subgraph, pos, edgelist=sigEdges, width=2.0, edge_cmap=corrMap,
                               edge_vmin=-1,
                               edge_vmax=1,                           
                               edge_color=[subgraph[e[0]][e[1]]['weight'] for e in sigEdges],
                               style='solid', alpha=1)    
    nx.draw_networkx_edges(subgraph, pos, edgelist=badEdges, width=1.0,                         
                               edge_color='grey',
                               style='solid', alpha=.3)    
    
    # add color bar as key to heats
    plt.colorbar()
    
    #nx.draw_networkx_edge_labels(subgraph,pos,edge_labels=eLab)
    if not args.nonames:
        nx.draw_networkx_labels(subgraph, pos, font_weight='bold', font_size=8)
    plt.axis('off')
    
    # write out the file(s)
    try:
        plt.savefig(args.out)
    except ValueError:
        plt.savefig('%s.png' % (args.out))
        
    if args.graphml:
        raise NotImplemented
        #nx.write_graphml(subgraph,args.graphml)
        
    if not args.noshow:
        plt.show()


if __name__ == "__main__":
    main()