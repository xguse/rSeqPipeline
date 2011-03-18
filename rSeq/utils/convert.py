import sys

def genephred2refflat(genePhredPath,refFlatPath):
	genePhred = open(genePhredPath,'rU')
	refFlat   = open(refFlatPath,'w')
	for line in genePhred:
		fields = line.rstrip('\n').split('\t')
		for i in range(len(fields)):
			fields[i] = fields[i].rstrip(',')
		if fields[0][-3] == "-":
			fields.insert(0,fields[0][:-3])
		else:
			fields.insert(0,fields[0])
		line = '\t'.join(fields)
		refFlat.write('%s\n' % (line))
	genePhred.close()
	refFlat.close()
		