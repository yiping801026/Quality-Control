import pysam
import numpy as np
import pandas as pd
import sys


def cal_GC(bam_file):
    sam = pysam.AlignmentFile(bam_file, 'rb')
    read_num=0
    GCnum=0
    base_num=0
    for read in sam:
        read_num+=1
        seq=read.query_sequence
        g_count=seq.count("G")
        c_count=seq.count("C")
        GCnum+=(g_count+c_count)
        #print('GCnum',GCnum)
        base_num +=len(seq)
        #print('base_num',base_num)
    GCcount = GCnum/base_num
    print("GC percentage of bam is ",GCcount)
    return GCcount


def main():
    cal_GC(sys.argv[1])

if __name__ == "__main__":
	main()

