from rSeq.utils.motifDiscovery.motifs import Motif,getHitDict
from gusPyCode.defs.bioDefs import makeRandomSeq
#m = Motif('WGATAR',mType='IUPAC')
#s = 'TGTACGAAGTGATAGTAGACTCGCTCCGCGAGCTTCTTACGTATAGCTCAAGGGTAAGCGGCAGTGTGGAAGTTCGATTTACGGGCTACGCGACCTGTTCCTACATCTCTGCAAGTATAATAAGATTCTGATTCCAGCAGTGATCCTGATCTATATGACGAATCAATGAGGTTAGGCAAGACTCATAGCGGTGGGGTAGTAAGCTTGTCATAACCAAGGACTGAATATTGCTCGGTCTGTATGCTACAATTTTTGGAACTTGCTTCCATACGTGCAACACGCACGTGAGGCTGGACGGGGATCCATGAAGGTGCGTGTCATGCAGCGCTTTGACGTGCCTTATCTCGAGCGTTAAGAGGATAAGCGAACACCAAATTGAGTCGGTTGAGCTATGGCCTCGAATCTCTACAGGCGCGGCCTAGTTACATTTGTTTCCTTCCGTAGTTCGCATAGGAAACTATGCGGGCAAGCAATGCAGCTCAACTCGGGCACATCGAGGTTCCCGGCAGCAAGTGATTGTCAGTGGTGTCAATGTCTCAGCAACCCCGTGGCAAGTCACTAGCCCAATGACAATTTCATTACACGAAATTGCTGATTCAGTCGTGGCCAGCCTCATTGCCTAGGTGCCTAATTCTCCAGGGTATAGGGTCAACTCGTCATAACGAAGCCATCGCTTACGTGGGGTATTCTGGTGCAAGCTGTAAGTTCGTCTCCTTTTCCTACACCCACCGTAGGGCAACTAAAGATAGATCGATCTGAGACTTTCTCGTTGCGAGTAGTCGGATTTTTATCAATTCATCCGCTTAGCGGCGGTTTCAGCGAGCATGCTGTGACCTAAAGGTGAGCTATGTGCATCCGCACCATTGGCTTTTTATACTTTAACTGGGGTGTGGGAACGCCACCCCCCAGGACTTCTCGCCTGGGGGAGCCAGAGTCGACGACCATACGTCTCCTGGATGGGATCCGTGTCGTGACATTGTATTTGGTAATCAAGCGCAGATGATTCCTCTGTCCAACCTCGCGTTTCGTGCTGTCACATGGGCTGGCGTGCGGGAAGACCCTGCTACACTCGGTGTATCTTTGTCAGTAACCATATGCGAACCTAGATCGAAATTACAGTACGCTGGCACAACACGAATCTCTACTATTCGAGGTTTGGGTGACCCCAATAGGATGCGCTAACAGGGGTAGGCTGCGCCCTCGGTAACTGGCCTAAACTGTCAAAGTTATCCGACATAGAGCCAGGCCTTAGTCATAACAGTAGATCCGCGGTACGGAAATCACACTCACATGACTCGATAAAGTGGGCGCAACAGCGGGCCATACCGTTGTGGCCCCCCGTTTCATGCCCCTTCGTTTGAGGTTGTGAGGACGAGGCTGTGCGTGGACAAGGCTGCAAACCGCCTACTTACGCCAGGGTCTGGCACTGCTAGAGGAACAATCAGGGGGTTACCGCTGTATTAAGGGATTCATTCCGTCCCACTGTCAGAAATCCCGTAGATAGCCCTAATGCTACATCGCTTTATGACTACTCATCTCCCGGCTGTCTACATAGGTCCCTGTGCATGGGCAGTGTCCTTTGGTGTCAGTGGTCCTCCCTGTATACCACGATTTGCGGTTTGCTCTAGGTGTACAGAAGTATTTATCTGTGTTGGTGGACCCTACGCCTGAATTTAGTCTGATACATCCTTGGTGCCCAAATCATTACATTATGTCGATATCACATAGTGGTGGTGAAGGGCTGGCAAAAGGTAAATCGGTATTGTCTTGCTGGGACAGCTCACTCCTCTGTTCGGAGACATTCCAGTTGCTGGAGCACTCTCAGTAAAGCCCGCTAAAGCGCCCGTTCGCATGTCCTAAGACGTCTGCACTCTATTCCGCTACCGGCAGCCCATGCGGGGTGTCCCTATCGTTGAAGACTCGGAAGGGATCACTCGCAATCAGACGCTGCATAAAAAGGATAAAAAGGGGTTTCTTCGTACCAGAATCACTCCCTGTCA'
#len(m.scan(s))
#len(m.scan(s))
#len(m.scan(s,pValThresh=0.1))
#len(m.scan(s,pValThresh=0.00000000000000000001))

d = {}
for i in range(1000):
    d[i] = makeRandomSeq(100)
  
mList = [Motif('WGATAR',mType='IUPAC'),Motif('WR',mType='IUPAC')]
hD = getHitDict(mList,d,pThresh=0.01,halfAT=0.25,halfGC=0.25)

None