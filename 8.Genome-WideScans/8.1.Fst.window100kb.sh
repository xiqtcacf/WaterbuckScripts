#!/bin/bash

batch=$1
export DIR=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc_poolpops_slidingwindow100k
cd $DIR

input1=/home/users/xi/Waterbuck_project/6.Population_structure_RM_6PK1K2_RM_AdmixedInds/poolingPop.saf.twoSubspecies/Goat_Chr
#file1=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc/$batch.1
#file2=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc/$batch.2
input2=/home/users/xi/Waterbuck_project/10.Fst/Goat_anc_poolpops_slidingwindow100k10k

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS

#for i in $file1;
# do
#	while read file ; do
#	   $angsd -bam $file.txt -out $file -doSaf 1 -anc $fasta -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2
#	   $realSFS $i $file -fold 1 -P 10 > $i.$file.folded.sfs
#	done < $test
#done

# Do 2dsfs from goat ancestral saf files, use it to estimate global fsts
while read file ; do

#$realSFS $input/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.$file.saf.idx $input/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.$file.saf.idx -P 30 > Common.Defassa.Goatanc.$file.2dsfs
#$realSFS fst index -whichFst 1 $input/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.$file.saf.idx $input/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.$file.saf.idx -P 30 -sfs Common.Defassa.Goatanc.$file.2dsfs -fstout Common.Defassa.Goatanc.$file
$realSFS fst stats2 $input2/Common.Defassa.Goatanc.$file.fst.idx -win 100000 -step 100000 > Common.Defassa.Goatanc.$file.fst.win100K

done < $input1/$batch
