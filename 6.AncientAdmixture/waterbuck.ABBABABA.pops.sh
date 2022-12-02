#!/bin/bash

file=$1
export DIR=/home/users/xi/Waterbuck_project/10.Abbababa/10pops_108inds.8outgroups/Defassa
cd $DIR

angsd=/home/users/xi/software/angsd/angsd
jackKnife=/home/users/xi/Waterbuck_project/10.Abbababa/jackKnife.R

fasta=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta
fai=/home/users/xi/Waterbuck_project/Reference/DefassaWaterbuck.fasta.fai
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

indname=/home/users/xi/Waterbuck_project/10.Abbababa/10pops_108inds.8outgroups/indNames.txt


$angsd -GL 2 -out Defassa.$file -nThreads 20 -doAbbababa 1 -doCounts 1 -bam $file -minMapQ 30 -minQ 20 -sites $sites_filtering -blockSize 5000000  -useLast 1
Rscript $jackKnife file=Defassa.$file.abbababa indNames=$indname outfile=Defassa.$file.abbababa boot=1

