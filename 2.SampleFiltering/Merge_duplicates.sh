#!/bin/bash

batch1=$1 ###merged bam name
batch2=$2 ###first bam need to be merged
batch3=$3 ###second bam need to be merged, and can be added more inds
batch4=$4 ###reference name 

export DIR=/home/users/xi/Waterbuck_project/7.merged_Samples
cd $DIR

input=/home/users/xi/Waterbuck_project/raw_bam/
output=/home/users/xi/Waterbuck_project/7.merged_Samples/
samtools=/home/users/xi/software/samtools-1.14/samtools
 
###### merge duplicated bams into one
$samtools merge $batch1.$batch4.bam $input/$batch2.$batch4.bam  $input/$batch3.$batch4.bam -@ 20

###### rename and index
$samtools reheader --command 'sed "s/\tSM:[0-9].../\tSM:Matetsi_1614_1614_1/g"' $batch1.$batch4.bam > $batch1.$batch4.Reheader.bam
$samtools index $output/$batch1.$batch4.Reheader.bam  

###### checking depth of merged samples
$samtools depth -a $batch1.$batch4.Reheader.bam > $batch1.$batch4.Reheader.bam.depth
cat $batch1.$batch4.Reheader.bam.depth | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $batch1.$batch4.Reheader.bam.depth.mappingCoverage
