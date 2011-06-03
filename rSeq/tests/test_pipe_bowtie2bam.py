from rSeq.utils.align import pipe_bowtie2bam

fastq1 = '/mnt/drobo01/solexa_xtra1b_rn/solexa_overflow/CBA/data4/Data/Intensities/BaseCalls/GERALD_04-12-2009_hnwitt/CBA_sequence.txt,/mnt/drobo01/solexa_xtra1b_rn/solexa_overflow/CBB/data3/Data/Intensities/BaseCalls/GERALD_08-12-2009_cnicolet/CBB_sequence.txt'


pipe_bowtie2bam('Aedes_aegypti.AaegL1.dna.toplevel',
                '/home/gus/tmp/test_pipe_bt2bam/test_CBx.bam',
                fastq1=fastq1,
                options='--solexa1.3-quals -v 2 -m 1 -S -p 2')