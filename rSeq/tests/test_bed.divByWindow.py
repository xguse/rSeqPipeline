from rSeq.utils.bed import divByWindow


bedA = '/Users/biggus/sandbox/testBEDtools/D.bed'
bedB = '/Users/biggus/sandbox/testBEDtools/C.bed'
outDir = '/Users/biggus/sandbox/testBEDtools/this/that'

result = divByWindow(bedA,bedB,win=[50,51],cols=[6,4],side='right',outDir=outDir)
