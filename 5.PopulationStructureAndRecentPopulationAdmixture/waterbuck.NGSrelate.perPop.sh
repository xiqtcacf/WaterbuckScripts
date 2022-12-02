#!/bin/bash

batch=$1
Hinds=$2
inds=$3
export DIR=/home/users/xi/Waterbuck_project/6.population_structure/NGSrelate/$batch
cd $DIR

angsd=/home/users/xi/software/angsd/angsd
ngsRelate=/home/users/xi/software/NGSrelate/ngsRelate/ngsRelate
bams=/home/users/xi/Waterbuck_project/6.population_structure/NGSrelate/$batch.txt
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
sampleID=/home/users/xi/Waterbuck_project/6.population_structure/NGSrelate/$batch.samplesID.list

### First we generate a file with allele frequencies (angsdput.mafs.gz) and a file with genotype likelihoods (angsdput.glf.gz).
$angsd -bam $bams -out waterbuck.$batch -GL 2 -sites $sites_filtering -doGlf 3 -minInd $Hinds -dosnpstat 1 -doHWE 1 -hwe_pval 1e-6 -doMajorMinor 1 -doMaf 1 -minmaf 0.05 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -nThreads 10
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat waterbuck.$batch.mafs.gz | cut -f5 | sed 1d > freq
### run NgsRelate
$ngsRelate -g waterbuck.$batch.glf.gz -f freq -n $inds -O waterbuck.$batch.res -z $sampleID -p 10
