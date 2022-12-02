#!/bin/bash

batch=$1
export DIR=/home/users/xi/Waterbuck_project/6.population_structure/SAF_eachInds_het_2DKing/DefassaWaterbuck/bam_list_perInd/
cd $DIR


angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
fasta=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta
fai=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta.fai
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions


while read file ; do

$angsd -bam $file.txt -out $file -doSaf 1 -anc $fasta -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2
$realSFS $file.saf.idx -P 3 -fold 1 > $file.saf.idx.folded.sfs

done < $batch.name 
