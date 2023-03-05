import pandas as pd
import sys


#### main
#### ------ read table ------ ####

flagstat = open(sys.argv[1]).readlines()

#### ------ get flagstat data ------ ####
alignment_num = int(flagstat[0].split(' ')[0])
paired_in_sequencing = int(flagstat[5].split(' ')[0])
PE = float(paired_in_sequencing)/float(alignment_num)
singletons = flagstat[10].split(' ')[4].replace("(","")
proper_pair = float(flagstat[8].split(' ')[0])/float(paired_in_sequencing)
duplicate_rate = float(flagstat[3].split(' ')[0])/alignment_num

print("Details of raw bam file",'\n',"alignment_num is",alignment_num,'\n',"PE is",PE,'\n',"proper_pair is ",proper_pair,"\n","duplicate_rate is",duplicate_rate)




