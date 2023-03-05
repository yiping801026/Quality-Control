### if loop for the input bed 4 types：
### 'all' for whole genome/gastric_colon for gastric&colon peaks/gastric/colon/gastric_filtered for selected peaks
### eg.: bash /home/pingyi/renkeQC_flag_plus/script/tss_peak.sh gastric ""

bed_choice=$1
input_bw_file=$2
output_file=$3
out_bed_file=$4

### 4 choices:
all=/home/pingyi/references/tss/refTSS_v3.3_human_coordinate.hg38.bed
gastric_colon=/home/pingyi/renkeQC_flag_plus/data/wo_data/all.txt
gastric=/home/pingyi/renkeQC_flag_plus/data/wo_data/gastric_probe.txt
colon=/home/pingyi/renkeQC_flag_plus/data/wo_data/colon_probe.txt
gastric_filtered=/home/pingyi/renkeQC_flag_plus/data/wo_data/filtered_gastric_probe.txt



### step 3. calculate the matrix && plot
if [ "$bed_choice" == "all" ];then
    computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -p 32 -R $all \
    -S $input_bw_file \
    --skipZeros -o $output_file --outFileSortedRegions $out_bed_file -q

elif [ "$bed_choice" == "gastric_colon" ];then
    computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -p 32 -R $gastric_colon \
    -S $input_bw_file \
    --skipZeros -o $output_file --outFileSortedRegions $out_bed_file -q



elif [ "$bed_choice" == "gastric" ];then
    computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -p 32 -R $gastric \
    -S $input_bw_file \
    --skipZeros -o $output_file --outFileSortedRegions $out_bed_file -q

elif [ "$bed_choice" == "colon" ];then
    computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -p 32 -R $colon \
    -S $input_bw_file \
    --skipZeros -o $output_file --outFileSortedRegions $out_bed_file -q
    
else
    computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -p 32 -R $gastric_filtered \
    -S $input_bw_file \
    --skipZeros -o $output_file --outFileSortedRegions $out_bed_file -q


fi