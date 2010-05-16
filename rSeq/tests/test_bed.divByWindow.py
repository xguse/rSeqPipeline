from rSeq.utils.bed import divByWindow


bedA = '/home/dunnw/data/genomes/aaegypti.SUPERCONTIGS-Liverpool.AaegL1/aaegypti.GENES-AaegL1.2.bed'
bedB = '/home/dunnw/data/solexa/bowtie_out/LSxLBx_bowtie.bed'
outDir = '/home/dunnw/data/solexa/bowtie_out/LSxLBx_bowtie_results'

print divByWindow(bedA,bedB,win=[2000,500],cols=[6,6],side='right',outDir=outDir)
