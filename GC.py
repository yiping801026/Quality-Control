import pysam
import numpy as np
import pandas as pd
import math
import sys
#import pyfaidx



def calGC(bamfile, bed):
	sampleid = bamfile.strip().split("/")[-1].split(".")[0]
	out = open("%s.GCcount.txt"%sampleid,"w")
	sam = pysam.AlignmentFile(bamfile, 'rb')
	bedfile = open(bed, "r")
	allgc = 0
	allbase = 0
	for lines in bedfile:
		chrs, start, end = lines.strip().split("\t")
		start = int(start)
		end = int(end)
		base_num = 0
		GCnum = 0
		read_num=0
		for read in sam.fetch(chrs,int(start),int(end)):
			if not read.is_duplicate and not read.is_secondary:
				read_num+=1
				seq = read.query_sequence
				seq_len = read.query_length
				read_start = read.reference_start +1
#				read_end = read.reference_end
				read_end = read_start + seq_len-1
				read_name = read.query_name
				if (read_start<= start) and ( start < read_end <= end):
					seq_in_bed = seq[start - read_start:]
				elif (read_start >= start) and (read_end <= end):
					seq_in_bed = seq
				elif ( start <= read_start <= end) and (read_end >= end):
					seq_in_bed = seq[:end-read_end-1]
				elif (read_start <= start) and (read_end >= end):
					seq_in_bed = seq[start-read_start:end-read_end-1]
				g_count = seq_in_bed.count("G")
				c_count = seq_in_bed.count("C")
				GCnum +=(g_count + c_count)
				base_num +=len(seq_in_bed)
				allgc +=GCnum
				allbase+=base_num
#		       	print(read_num,GCnum,base_num)
		if base_num==0:
			GCcount=0
		else:
			GCcount = GCnum / base_num
		out.write(lines.strip() +"\t"+str(GCcount)+"\t"+str(read_num)+"\n")
	gchan = allgc/allbase*100
	out2 = open("%s.all_GCcount.txt"%sampleid,"w")
	out2.write(sampleid +"\t"+str(gchan)+"\n")
	print('GC_wotss is',gchan)
#	return sampleid,dict1
	return

#calGC("/home/pingyi/renkeQC_flag_plus/data/wo_data/ZP-39_BTSB0P0006-1-MPN1_sorted.bam","/home/pingyi/renkeQC_rm_plot/data/sorted.bed")


def main():
	### 
	calGC(sys.argv[1],sys.argv[2])


if __name__ == "__main__":
	main()
