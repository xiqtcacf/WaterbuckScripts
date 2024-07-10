#!/bin/bash

#batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_NullModel_NonParametricBootstraps
cd $DIR

fasta=/home/users/xi/Waterbuck_project/10.fastsimcoal27/Red_lechwe.zoo100.DefassaWaterbuck.fsc27.fa.gz
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/DefassaWaterbuck.GT100k.sites_filter/DefassaWaterbuck.GT100k.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

Goat_fasta=/home/users/xi/Waterbuck_project/10.fastsimcoal27/Red_lechwe.zoo100.Goat.fsc27.fa.gz
Goat_sites_filtering=/home/users/xi/Waterbuck_project/10.sex_chromosome_goat/Goat.allChr.all_rep_het_dep_map.regions

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
winsfs=/davidData/malthe/projects/cov/sim/bin/winsfs

### get 100 bootstraps 2dsfs using a genomic jackknife based on 100 blocks by winsfs
$winsfs split --threads 40 -S 100 --tolerance 1e-8 --method pseudo-loo --sfs ../WinsfsNeed.Matetsi.Samole.Red_lechweANC.DefassaWaterbuck.2dsfs ../Matetsi.Red_lechwe.DefassaWaterbuck.saf.idx ../Samole.Red_lechwe.DefassaWaterbuck.saf.idx  > winSFS.Matetsi.Samole.Red_lechweANC.DefassaWaterbuck.2dsfs.100bootstraps

### transfer winsfs 2dsfs 100bootstraps into fsc27 input file, first need to split per line as individual file
input=/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_NullModel_NonParametricBootstraps/get_InputFromWinsfs
for i in {100..200}
do
Rscript $input/get.fsc27Input.From2Dsfswinsfs.R $input $i
cat $input/header.row $input/bootstrap.2dsfs.matrix > $input/y
paste $input/header.column $input/y > $input/$i.2dsfs 
rm $input/bootstrap.2dsfs.matrix $input/y
done
