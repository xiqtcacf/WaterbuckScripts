#!/bin/bash

export DIR=/home/users/xi/Waterbuck_project/10.treemix/IBD_NJtree
cd $DIR

#pcangsd=/home/users/xi/software/pcangsd/pcangsd.py
bam=/home/users/xi/Waterbuck_project/10.treemix/IBD_NJtree/108inds.8outgroups.singleReadSampling.For_treemix.bam_full_path.txt
angsd=/home/users/xi/software/angsd/angsd
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions


###### calculate observation txt
#for i in {1..18}
#do
#   while read file ; do
#       eval vcftools --bcf $bcf/chr${i}.bcf.gz --indv $file --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
#       python $scripts/FiltervariantsNOancestral.py $freq/outgroup_chr${i}.freqs 1000 $txt/weights_chr${i}.txt Xi_input_observations/$file.chr${i}.observations_NOancestral.txt
#       done < $ind
#done

### single read sampling with angsd
$angsd -bam $bam -sites $sites_filtering -out 108inds.8outgroups.singleReadSampling -minMapQ 30 -minQ 20 -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -minMaf 0.05 -P 20


