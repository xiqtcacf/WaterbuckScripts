#!/bin/bash

#batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.Abbababa2/
cd $DIR



angsd=/home/users/xi/software/angsd/angsd
fasta=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta
fai=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta.fai
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
input=/home/users/xi/Waterbuck_project/outgroup/

while read file ; do

  $angsd -doFasta 2 -doCounts 1 -i $input/$file -minQ 20 -minMapQ 30 -P 20 -out $file.abbababa -sites $sites_filtering

done < outgroups.bams.txt
