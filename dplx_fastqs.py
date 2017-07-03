import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline

from dms2dfe.lib.io_seq_files import *

def contains_barcode(barcode_seq, read_seq, pos='start') : 
    """
    This searches for barcode sequence in given sequence read \
    in windows of +/- 5 nucleotide upstream and downstream.
    
    :param barcode_seq: sequence of barcode (string).
    :param barcode_seq: target sequence (string).
    """
    for i in range(0, 5) : # (12,17) : #start looking around the end of the 15th base (4 char tag + 11 char mid)
        if pos=='start':
            if barcode_seq in str(read_seq[i:]) :
                return True
        if pos=='end':
            if barcode_seq in str(read_seq[:i]) :
                return True
    if pos=='any':
#         print '%s\n%s' % (read_seq,barcode_seq)
        if barcode_seq in read_seq:
            return True

def R12merged(fastq_R1_read,fastq_merged_f):
    """
    This merges and joins R1 and R2 reads.
    
    merged files (*_dplxd_merged.fastq) have the R1 and R2 mixed together without joining. 
    eg. R1 R2 R1 R2 ..
    merged files (*_dplxd_joined.fastq) have the R1 and R2 joined together.
    eg.	R1R2 R1R2 ..
            
    :param fastq_R1_read: sequence of R1 read of .fastq (str).
    :param fastq_R2_read: sequence of R2 read of .fastq (str).
    :param fastq_merged_f: file object for output merged .fastq.  
    :param fastq_joined_f: file object for output joined .fastq.
    """
    from dms2dfe.lib.convert_seq import revcom,revers
    #MERGED
    fastq_R1_seq=str(fastq_R1_read.seq)
    fastq_R1_qua=''.join([chr(x + 33) for x in fastq_R1_read.letter_annotations['phred_quality']])
    fastq_writer(fastq_merged_f,fastq_R1_read.id+"R1",fastq_R1_seq,fastq_R1_qua)

def fastq2dplx_one_sided(fastq_R1_fh,barcode_R1s,fastq_fns,force=False) :
    """
    This demultiplexes fastq.
    
    :param fastq_R1_fh: path to R1 .fastq file.
    :param barcode_R1s: barcodes to demultiplex forward reads.
    :param fastq_fns: file names of demultiplexed .fastq files (saved in same folder as the non-multiplexed input file).
    """
    from dms2dfe.lib.convert_seq import cmplmt,revers,revcom
    
    bars2fns=dict(zip(barcode_R1s,fastq_fns))
    
    fastq_R1_reads=SeqIO.parse(fastq_R1_fh, 'fastq')
    fastq_dh=dirname(fastq_R1_fh)
    for fastq_R1_read in fastq_R1_reads :
        for barcode_R1 in bars2fns: #zip(barcode_R1s,barcode_R2s,fastq_fns):
#             print barcode_R1
            if contains_barcode(barcode_R1,
                                str(fastq_R1_read.seq) ,
                               pos='any'):
                fastq_merged_fh="%s/%s_dplxd_merged.qcd.fastq" % (fastq_dh,bars2fns[barcode_R1])
                with open(fastq_merged_fh,'a') as fastq_merged_f:
#                     print fastq_merged_fh
                    R12merged(fastq_R1_read,fastq_merged_f)
    fastq_R1_reads.close()

# fastq2dplx_one_sided(fastq_R1_fh,barcode_R1s,fastq_fns,force=True)