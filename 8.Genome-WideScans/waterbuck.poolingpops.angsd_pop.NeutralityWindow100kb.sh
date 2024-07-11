#!/bin/bash

batch=$1
Goat_Chr=$2
export DIR=/home/users/xi/Waterbuck_project/10.Neutrality_test_poolingGoat/
cd $DIR

input=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/poolingPop.saf.twoSubspecies/Goat_Chr

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
thetaStat=/home/users/xi/software/angsd/misc/thetaStat

while read file ; do

$realSFS $input/$batch.$file.saf.idx -P 30 > $batch.$file.saf.idx.sfs

$realSFS saf2theta $input/$batch.$file.saf.idx -outname $batch.$file -sfs $batch.$file.saf.idx.sfs
$thetaStat do_stat $batch.$file.thetas.idx -win 100000 -step 100000  -outnames $batch.$file.thetas.win100k.thetasWindow.gz

done < $input/$Goat_Chr
