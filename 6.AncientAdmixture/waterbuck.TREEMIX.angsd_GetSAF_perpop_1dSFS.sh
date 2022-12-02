#!/bin/bash

#batch=$1
export DIR=/brendaData/xi/African1Kg/Waterbuck/6.Population_structure_RM_6PK1K2_RM_AdmixedInds
cd $DIR


angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
fasta=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta
fai=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta.fai
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
input=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/pop.108inds.no_outgroup_crashed_duplicate_within_samples_Waterbuck.nok1k2_5inds.noAdmix12inds


while read file ; do

$angsd -bam $input/$file.list -out $file.Defanc -doSaf 1 -anc $fasta -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2 -nThreads 30
$realSFS $file.Defanc.saf.idx -P 30 > $file.Defanc.saf.idx.sfs

done < $input/pop.list
