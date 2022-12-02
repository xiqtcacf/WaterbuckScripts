#!/bin/bash

batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc
cd $DIR

input=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc/input
file1=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc/$batch.1
file2=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc/$batch.2

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS

# Do 2dsfs from goat ancestral saf files, use it to estimate global fsts
while IFS= read -r line1 && IFS= read -r line2 <&3; do

   $realSFS $input/$line1 $input/$line2 -P 30 > $line1.$line2.Goatanc.2dsfs
   $realSFS fst index -whichFst 1 $input/$line1 $input/$line2 -P 30 -sfs $line1.$line2.Goatanc.2dsfs -fstout $line1.$line2.Goatanc
   $realSFS fst stats $line1.$line2.Goatanc.fst.idx > $line1.$line2.Goatanc.Fst

done < $file1 3< $file2
