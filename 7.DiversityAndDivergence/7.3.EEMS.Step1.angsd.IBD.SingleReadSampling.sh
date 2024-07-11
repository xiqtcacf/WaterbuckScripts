#!/bin/bash

export DIR=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/SingleReadSampling_EEMS_input/
cd $DIR

#pcangsd=/home/users/xi/software/pcangsd/pcangsd.py
bam=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/SingleReadSampling_EEMS_input/108inds.no_outgroup_crashed_duplicate_within_samples_Waterbuck.nok1k2_5inds.noAdmix12inds.raw_bam_full_path.txt
angsd=/home/users/xi/software/angsd/angsd
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

### single read sampling with angsd
$angsd -bam $bam -sites $sites_filtering -out RM_relate5Inds.RMadmix12inds.waterbuck.DefassaWaterbuckRef.SingleReadSampling.ForEEMS -minMapQ 30 -minQ 20 -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.05 -P 30
