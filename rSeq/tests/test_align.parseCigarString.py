from rSeq.utils.align import cigarStr2AlignCoords
cig = '2MD3M6I2D2M'
print cigarStr2AlignCoords(cig,0,13,'+',intron='I')


