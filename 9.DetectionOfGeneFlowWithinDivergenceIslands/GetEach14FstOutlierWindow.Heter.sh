#!/bin/bash

#batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/add_NGSadmix3Inds_14region_heter
cd $DIR

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS

bam=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/info/nostar.add3NGSadmixInds.GoatBam.111inds.no_outgroup_crashed_duplicate_within_samples_Waterbuck.nok1k2_5inds.noAdmix12inds.raw_bam_full_path.txt
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
bed=/home/users/xi/Waterbuck_project/10.NJtree_PCAngsd_FstTwisstSlidingWindow/info/FstWindow14.bed
input=/home/mikkel/projects/2103_africa1k/waterbuck/filtered_bams/output/

fasta=/home/users/xi/Waterbuck_project/Reference/Goat.fasta
fai=/home/users/xi/Waterbuck_project/Reference/Goat.fasta.fai

while read chr ; do
 while read file ; do
  $angsd -i $input/$file -sites $sites_filtering -r $chr -out $chr.$file -anc $fasta -doSaf 1 -minMapQ 30 -minQ 20 -GL 2 -nThreads 20
  $realSFS $chr.$file.saf.idx -P 20 > $chr.$file.saf.idx.sfs
  done < $bam
done < $bed
