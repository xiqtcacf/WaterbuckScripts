#!/bin/bash

#### bash waterbuck.poolingpops.angsd_ngsLD.sh common_bula_bula Common 41 0.01

name=$1
species=$2
n_ind=$3
rnd_sample=$4
export DIR=/home/users/xi/Waterbuck_project/10.ngsLD/$species
cd $DIR

Goat_fasta=/home/users/xi/Waterbuck_project/Reference/Goat.fasta
Goat_fai=/home/users/xi/Waterbuck_project/Reference/Goat.fasta.fai
Goat_sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
Goat_Chr=/home/users/xi/Waterbuck_project/10.ngsLD/Chr.onlyFstOutliers.list

input=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/pop.108inds.no_outgroup_crashed_duplicate_within_samples_Waterbuck.nok1k2_5inds.noAdmix12inds
angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
thetaStat=/home/users/xi/software/angsd/misc/thetaStat
ngsLD=/home/users/xi/software/ngsLD/ngsLD

while read file ; do
$angsd -bam $input/$name.Goat.list -ref $Goat_fasta -out $name.GoatRef.$file -r $file -sites $Goat_sites_filtering  \
        	-minMapQ 30 -minQ 20 -doCounts 1 \
        	-GL 2 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        	-doGlf 4 -SNP_pval 1e-6 -nThreads 30

zless $name.GoatRef.$file.mafs.gz | cut -f1,2 > $name.GoatRef.$file.mafs.pos ####get pos file
NSITES=`zcat $name.GoatRef.$file.mafs.gz | tail -n+2 | wc -l` ### how to calculate the number of sites
$ngsLD --geno $name.GoatRef.$file.glf.gz --min_maf 0.05 \
       --posH $name.GoatRef.$file.mafs.pos --extend_out --log_scale T --n_threads 30 --n_ind $n_ind --n_sites $NSITES --rnd_sample $rnd_sample \
       --outH maf0.05.$name.GoatRef.$file.$rnd_sample.$n_ind.ld
done < $Goat_Chr
