#!/bin/bash

export DIR=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/no_NGSadmixInds_14region_PCAngsd
cd $DIR

pcangsd=/home/users/xi/software/pcangsd/pcangsd.py
angsd=/home/users/xi/software/angsd/angsd

bam=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/info/GoatBam.108inds.no_outgroup_crashed_duplicate_within_samples_Waterbuck.nok1k2_5inds.noAdmix12inds.raw_bam_full_path.txt
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
bed=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/info/FstWindow14.bed

while read file ; do
#### generate genotype likelihood of Beagles
$angsd -bam $bam -sites $sites_filtering -r $file -out $file.RM_relate5Inds.RMadmix12inds.waterbuck.GoatRef.ForPCAngsd -minMapQ 30 -minQ 20 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05  -nThreads 20
#### do PCAangsd
 for i in {1..10}
  do
  python3 $pcangsd -beagle $file.RM_relate5Inds.RMadmix12inds.waterbuck.GoatRef.ForPCAngsd.beagle.gz -e $i -o $file.RM_relate5Inds.RMadmix12inds.waterbuck.GoatRef.ForPCAngsd.Pcangsd.e${i} -minMaf 0.05 -threads 20
  done
done < $bed

