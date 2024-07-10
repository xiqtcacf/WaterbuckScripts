#!/bin/bash

batch=$1
export DIR=/home/wlk579/0.waterbuck/Fst_WaterC
cd $DIR

file1=/home/wlk579/0.waterbuck/Fst_WaterC/$batch.1
file2=/home/wlk579/0.waterbuck/Fst_WaterC/$batch.2

angsd=/home/wlk579/Server_bos/apps/angsd/angsd
realSFS=/home/wlk579/Server_bos/apps/angsd/misc/realSFS
input=/home/wlk579/0.waterbuck/Fst_WaterC/Defanc.saf.perPop

# Do 2dsfs from waterbuck saf files, use it to estimate global fsts
while IFS= read -r line1 && IFS= read -r line2 <&3; do

   $realSFS $input/$line1 $input/$line2 -P 20 > $line1.$line2.Defanc.2dsfs
   $realSFS fst index -whichFst 1 $input/$line1 $input/$line2 -P 20 -sfs $line1.$line2.Defanc.2dsfs -fstout $line1.$line2.Defanc
   $realSFS fst stats $line1.$line2.Defanc.fst.idx > $line1.$line2.Defanc.Fst

done < $file1 3< $file2
