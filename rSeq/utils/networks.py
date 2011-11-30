import os,sys

try:
    from xml.etree import cElementTree
except AttributeError:
    pass
import networkx as nx

import scipy


from rSeq.utils.files import tableFile2namedTuple


def import_edges(graphObj,edgeTablePath,startNodeHeader,endNodeHeader):
    """
    GIVEN:
    1) graphObj = networkx.Graph() instance
    2) edgeTablePath = path to labeled tsv file containing the edge relationships.
    3) startNodeHeader = header label containing the first node of the edge.
    4) endNodeHeader = header label containing the second node of the edge.


    DO:
    1) Populate graphObj with the edges.
    
    RETURN:
    1) None
    """
    
    edgeTable = tableFile2namedTuple(edgeTablePath)

    graphObj.add_edges_from([(edge.__getattribute__(startNodeHeader),edge.__getattribute__(endNodeHeader)) for edge in edgeTable])
    
    
def graph_connected_nodes(graphObj,nodeList):
    """
    GIVEN:
    1) graphObj = networkx.Graph() instance.
    2) nodeList = list of nodes for which connected nodes should be extracted.

    DO:
    1) Extract any node in graphObj that shares and edge with any node in nodeList.
    2) Create subgraph described by the extracted nodes.
    
    RETURN:
    1) subgraph
    """
    
    # hash nodeList with set() for faster searching and query the set with
    # nodes in edges.
    nodes = set(nodeList)
    subNodes = []
    for node1,node2 in graphObj.edges_iter():
        if (node1 in nodes) or (node2 in nodes):
            subNodes.extend([node1,node2])
            
    return nx.Graph(graphObj.subgraph(subNodes))

def weight_edges_with_pearsonr(graphObj,dataVectors,uni=False):
    """
    GIVEN:
    1) graphObj = networkx.Graph() instance.
    2) dataVectors = dict of expression vectors with keys==TxName
    3) uni = if true, weight all edges the same
    
    DO:
    1) Iterate through all edges in graphObj.
    2) Calculate pearsonr bt nodes and convert this to length of edge.
    3) Overwrite edge and add 'length' attrib.
    
    RETURN:
    1) None
    """
    if not uni:
        for node1,node2 in graphObj.edges_iter():
            try:
                r_val,p_val = scipy.stats.pearsonr(dataVectors[node1],dataVectors[node2])
                graphObj.add_edge(node1, node2, attr_dict={'rVal':r_val,'pVal':p_val,'weight':r_val})
            except KeyError as e:
                sys.stderr.write('node not found: %s\n' % str(e))
    else:
        for node1,node2 in graphObj.edges_iter():
            try:
                r_val,p_val = scipy.stats.pearsonr(dataVectors[node1],dataVectors[node2])
                graphObj.add_edge(node1, node2, attr_dict={'rVal':r_val,'pVal':p_val,'weight':0})        
            except KeyError as e:
                sys.stderr.write('node not found: %s\n' % str(e))            