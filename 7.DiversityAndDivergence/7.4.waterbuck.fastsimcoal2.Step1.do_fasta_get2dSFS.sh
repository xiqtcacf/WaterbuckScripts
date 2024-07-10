#!/bin/bash

#batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.fastsimcoal27
cd $DIR

input=/home/users/xi/Waterbuck_project/10.sex_chromosome_goat/$batch.list

fasta=/home/users/xi/Waterbuck_project/10.fastsimcoal27/Red_lechwe.zoo100.DefassaWaterbuck.fsc27.fa.gz
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

Goat_fasta=/home/users/xi/Waterbuck_project/10.fastsimcoal27/Red_lechwe.zoo100.Goat.fsc27.fa.gz
Goat_sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS

### do_fasta, using the red lechwe (Kobus leche) to represent the ancestral allele.
while read file ; do
$angsd -doFasta 2 -doCounts 1 -i $input/Red_lechwe.zoo100.DefassaWaterbuck.bam -minQ 20 -minMapQ 30 -P 20 -out Red_lechwe.zoo100.DefassaWaterbuck.fsc27 -sites $sites_filtering
done < outgroups.bams.txt

### get saf per pop and the globle 2dsfs
$angsd -bam Matetsi.list -out Matetsi.Red_lechwe.DefassaWaterbuck -doSaf 1 -anc $fasta -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2 -nThreads 30
$angsd -bam Samole.list -out Samole.Red_lechwe.DefassaWaterbuck -doSaf 1 -anc $fasta -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2 -nThreads 30
$realSFS Matetsi.Red_lechwe.DefassaWaterbuck.saf.idx Samole.Red_lechwe.DefassaWaterbuck.saf.idx  -P 30 > Matetsi.Samole.Red_lechweANC.DefassaWaterbuck.2dsfs
