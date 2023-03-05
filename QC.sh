################## input files #####################
#script
#GC_wotss_script=/home/pingyi/renkeQC_rm_plot/script/GC.py 
#tss_peak_script=/home/pingyi/renkeQC_flag_plus/script/tss_peak.sh

#for all
tss_bed="/home/pingyi/renkeQC_rm_plot/data/sorted.bed"
bed_path=/home/pingyi/renkeQC_flag_plus/script/NanOnCT_Panel_v1.0_Covered_hg38.bed
bed_file=/home/pingyi/references/tss/hg38-blacklist.v2.bed

#for indiuvidual
sample_name=XHRK0P0002-3-MPN1
dirin_file=/home/pingyi/sample_pool/XHRK0P0002-3-MPN1
bam_file="$dirin_file/XHRK0P0002-3-MPN1.bam"
sort_bam_file="$dirin_file/XHRK0P0002-3-MPN1_sorted.bam"
### output files
#dirout_file=/home/pingyi/renkeQC_flag_plus/out/wo_out/wo_out_test
dirout_file=/home/pingyi/test

################## for '.fastq' #####################

### 1.base_sum || 7. len_avg
cd $dirin_file
seqkit stat $dirin_file/"${sample_name}_1.fq.gz" $dirin_file/"${sample_name}_2.fq.gz"
### 15. GC: 
unzip $dirin_file/"${sample_name}_1_fastqc.zip" -d $dirin_file/"tmp"
GC_1=`awk 'NR==10{print $2}' $dirin_file/tmp/"${sample_name}_1_fastqc"/fastqc_data.txt`
unzip $dirin_file/"${sample_name}_2_fastqc.zip" -d $dirin_file/"tmp"
GC_2=`awk 'NR==10{print $2}' $dirin_file/tmp/"${sample_name}_2_fastqc"/fastqc_data.txt`
echo "raw GC of "${sample_name}_1.fq.gz" is $GC_1"
echo "raw GC of "${sample_name}_2.fq.gz" is $GC_1"
rm -r $dirin_file/"tmp"

################## for '.bam' #####################

### 7. length_ave
samtools view $bam_file | \
awk 'BEGIN{length_sum;a}{a+=1;length_sum+=length($10)}END{print "length_ave is" "\t" length_sum/a}' 
wait
### index for sorted.bam 
samtools index -b $sort_bam_file
wait
####### calculate depth for sorted.bam of every regions
samtools depth -aa $sort_bam_file > $dirout_file/"${sample_name}.sorted.depth"  
### 8. depth 
depth=`awk 'BEGIN{sum;a;b}{if($3>0)a+=1;sum+=1;b+=$3}END{print "depth_ave is""\t"b/sum}' $dirout_file/"${sample_name}.sorted.depth" `
echo "depth is : $depth"
### 9. coverage 
coverage=`awk 'BEGIN{sum;a;b}{if($3>0)a+=1;sum+=1;b+=$3}END{print "coverage is""\t"a/sum}' $dirout_file/"${sample_name}.sorted.depth" ` 
echo "coverage is : $coverage"
### 10. coverage_6x
coverage_6x=`awk 'BEGIN{sum;a;b}{if($3>6)a+=1;sum+=1;b+=$3}END{print "coverage_6x is""\t"a/sum}' $dirout_file/"${sample_name}.sorted.depth"`
echo "coverage_6x is : $coverage_6x"
### 11. coverage_3x
coverage_3x=`awk 'BEGIN{sum;a;b}{if($3>3)a+=1;sum+=1;b+=$3}END{print "coverage_3x is""\t"a/sum}' $dirout_file/"${sample_name}.sorted.depth"`
echo "coverage_3x is : $coverage_3x"
### 17. tss_peak_plot
bedtools intersect -ubam -a $sort_bam_file -b $bed_file -v > $dirout_file/"${sample_name}.filtered.bam"  
wait
# step 1. add index to the bam files 
samtools index $dirout_file/"${sample_name}.filtered.bam"  
wait
#step 2. create the .bw file from .bam files
bamCoverage -b $dirout_file/"${sample_name}.filtered.bam"  -o $dirout_file/"${sample_name}.filtered.bw" --minFragmentLength 120 --maxFragmentLength 200
#echo "Done for step3" 
wait
### tss_peak_plot: have 5 types: all/gastric_colon/gastric/colon/gastric_filtered,eg.: bash /home/pingyi/renkeQC_flag_plus/script/tss_peak.sh gastric
bash /home/pingyi/renkeQC_flag_plus/script/tss_peak.sh all $dirout_file/"${sample_name}.filtered.bw" $dirout_file/"${sample_name}.2k_TSS.gz" $dirout_file/"${sample_name}.2k_TSS.out.bed"
#echo "Done for step4" 
wait
plotHeatmap -m $dirout_file/"${sample_name}.2k_TSS.gz" -o $dirout_file/"${sample_name}.2k_TSS.png"

### get the 4.alignment reads 5.PE 6. proper_pair 12. duplicate_rate
samtools flagstat $bam_file > $dirout_file/"${sample_name}_flagstat.txt"
python /home/pingyi/renkeQC_flag_plus/script/flagread.py $dirout_file/"${sample_name}_flagstat.txt"

### 15. GC GC% of bam file 
python /home/pingyi/renkeQC_rm_plot/script/rm_tss_bed.py $sort_bam_file

### 16. GC_wotss: GC% without tss regions 
python /home/pingyi/renkeQC_rm_plot/script/GC.py $sort_bam_file $tss_bed


