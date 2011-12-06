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
import scipy.stats as stats

from rSeq.utils.errors import *
from rSeq.utils.files import tableFile2namedTuple
from rSeq.utils.networks import import_edges
from rSeq.utils.networks import graph_connected_nodes
from rSeq.utils.networks import weight_edges_with_pearsonr
from rSeq.utils.expression import mangle_expn_vectors

from rSeq.scripts.pearson_scatter import plotScatter



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
    
    parser.add_argument('--cmap', type=str, required=True, nargs='+',
                        help="""A list of species-prefix:color combinations to set the node colors:
[--cmap <AAEL:b ...>]. the number of combinations should match the number of files given to --xprnfile.""")
    
    parser.add_argument('--log', action='store_true',
                        help="""Plot the points on a log:log scale. (Default: %(default)s)""")
    parser.add_argument('--show', action='store_true',
                        help="""Plot the image for interactive manipulation, otherwise just write the file. (Default: %(default)s)""")
    parser.add_argument('--pdf', action='store_true',
                        help="""Plot the image as a pdf: png otherwise. Png is preferable when data size is large. (Default: %(default)s)""")
    parser.add_argument('--out', type=str, default='',
                        help="""Base path for output. (Default: current working directory)""")
    parser.add_argument('--pdf', type=str, default=False,
                        help="""Load graph from a gpickle. (Default: %(default)s)""")
    
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
    
    load_pickle = True
    if load_pickle:
        subgraphs = nx.read_gpickle('/tmp/ortho_weighted_subgraphs.gpickle')
    else:
        # lets get started: init the graph
        graph = nx.Graph()
        
        for f in args.hfile:
            import_edges(graphObj=graph,edgeTablePath=f[0],startNodeHeader=f[1],endNodeHeader=f[2])
        
            
        # remove the '' node caused by unpaired relationships
        graph.remove_node('')
        
    
        # weight the edges in each graph by the pearsonr between their expression vectors
        weight_edges_with_pearsonr(graphObj=graph,dataVectors=xDict,uni=False)
        
            
        # if the edge length is imposible to graph (inf or nan) kill the edge
        #badEdges = []
        #edgesMissingNodes = []
        #for i,j in graph.edges_iter():
            #try:
                #if math.isnan(graph[i][j]['rVal']) or math.isinf(graph[i][j]['rVal']):
                    #badEdges.append((i,j))
            #except KeyError:
                #edgesMissingNodes.append((i,j))
                    
        #graph.remove_edges_from(badEdges)
        #graph.remove_edges_from(edgesMissingNodes)
        
        # Get all subgraphs
        subgraphs = nx.connected_component_subgraphs(graph)  
        nx.write_gpickle(subgraphs,"/tmp/ortho_weighted_subgraphs.gpickle")
        print "I layed a pickle!!"
    
    args.galaxy = False
    #args.label2 = "Pct w/ significant positive corr (r >= 0.5, p <= 0.05)"
    args.label2 = "Pct with significant correlation" 
    for prefix in cDict:
        args.label1 = "Usable paralogs per subgraph within %s" % (prefix)
        #args.label1 = "%s x" % (prefix) 
        pearsonStats,data = get_within_data(prefix,subgraphs)
        plotScatter(pearsonStats,data,args,color=cDict[prefix])
        
    args.label1 = "Usable orthologs per subgraph between AGAP and CPIJ" 
    #args.label1 = "both x"  
    pearsonStats,data = get_between_data(prefixes=cDict.keys(),subgraphs=subgraphs)
    plotScatter(pearsonStats,data,args,color='green')
    
    print "Done."
    

def get_within_data(prefix,subgraphs):
    """
    """
    x_data = []
    y_data = []
    
    for sGraph in subgraphs:
        usableEdges = 0.0
        goodEdges = 0.0
        for node1,node2 in sGraph.edges_iter():
            if node1.startswith(prefix) and node2.startswith(prefix):
                usableEdges += 1
                try:
                    if (sGraph[node1][node2]['rVal'] >= 0.8) and (sGraph[node1][node2]['pVal'] <= 0.05):
                        goodEdges += 1
                except KeyError:
                    usableEdges -= 1
        
        try:
            y_data.append((goodEdges/usableEdges)*100)
            x_data.append(usableEdges)
        except ZeroDivisionError:
            #print "good:%s usable:%s" % (goodEdges,usableEdges)
            pass
    
    x_data,y_data = tuple(x_data),tuple(y_data)
    pStats = stats.pearsonr(x_data,y_data)
    
    return (pStats,(x_data,y_data))
    
def get_between_data(prefixes,subgraphs):
    """
    """
    x_data = []
    y_data = []
    
    for sGraph in subgraphs:
        usableEdges = 0.0
        goodEdges = 0.0
        for node1,node2 in sGraph.edges_iter():
            testStr = str(node1)+str(node2)
            if (prefixes[0] in testStr) and (prefixes[1] in testStr):
                usableEdges += 1
                try:
                    if (sGraph[node1][node2]['rVal'] >= 0.8) and (sGraph[node1][node2]['pVal'] <= 0.05):
                        goodEdges += 1
                except KeyError:
                    usableEdges -= 1
        
        try:
            y_data.append((goodEdges/usableEdges)*100)
            x_data.append(usableEdges)
        except ZeroDivisionError:
            #print "good:%s usable:%s" % (goodEdges,usableEdges)
            pass
    
    x_data,y_data = tuple(x_data),tuple(y_data)
    pStats = stats.pearsonr(x_data,y_data)
    
    return (pStats,(x_data,y_data))    

if __name__ == "__main__":
    main()