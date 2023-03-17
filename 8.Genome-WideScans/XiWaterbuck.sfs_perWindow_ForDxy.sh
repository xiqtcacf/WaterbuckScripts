#! /bin/bash

r_chr_bed=$1
chr=$2

#####bash XiWaterbuck.sfs_perWindow_ForDxy.sh NC_030808.1.rFormat.100Kwindow.Goat.fasta.chr29.bed chr1

export DIR=/home/users/xi/Waterbuck_project/0_add_Dxy/asfsp
cd $DIR

angsd=/home/users/xi/software/angsd/angsd
realSFS=/home/users/xi/software/angsd/misc/realSFS
BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools

output=/home/users/xi/Waterbuck_project/0_add_Dxy/dxy_allAutoChr/$chr

fasta=/brendaData/xi/African1Kg/Waterbuck/Reference/Goat.fasta
sites_filtering=/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
Common_bam=/home/users/xi/Waterbuck_project/0_add_Dxy/pop_chr_info/dxy.Common.list
Defassa_bam=/home/users/xi/Waterbuck_project/0_add_Dxy/pop_chr_info/dxy.Defassa.list

#### step1 get window size bed for each chromosome and change format from bed to angsd input '-r chr2:135000000-140000000'
$BEDTOOLS makewindows -g Goat.fasta.chr29.txt -w 100000 > 100Kwindow.Goat.fasta.chr29.bed ###sutosome
awk '{print $1":"$2"-"$3}' 100Kwindow.Goat.fasta.chr29.bed > rFormat.100Kwindow.Goat.fasta.chr29.bed
$BEDTOOLS makewindows -g Goat.fasta.chrSex.txt -w 100000 > 100Kwindow.Goat.fasta.chrSex.bed ####sex chromosome
awk '{print $1":"$2"-"$3}' 100Kwindow.Goat.fasta.chrSex.bed > rFormat.100Kwindow.Goat.fasta.chrSex.bed

#### step2 using angsd to calculate saf per subspecies per window
while read window_bed; do
      $angsd -bam $Common_bam -out Common.$window_bed -doSaf 1 -anc $fasta -r $window_bed -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2
      $angsd -bam $Defassa_bam -out Defassa.$window_bed -doSaf 1 -anc $fasta -r $window_bed -sites $sites_filtering -minMapQ 30 -minQ 20 -GL 2
      $realSFS Common.$window_bed.saf.idx Defassa.$window_bed.saf.idx -P 20 -fold 1 > Common.Defassa.$window_bed.2dsfs
      rm Common.$window_bed.saf.*
      rm Defassa.$window_bed.saf.*
      rm *.arg

#### step3 calculate dxy per window
      python3 asfsp.py -I Common.Defassa.$window_bed.2dsfs -D 82,134 --calc dxy > $window_bed.dxy
      mv Common.Defassa.$window_bed.2dsfs ../dxy_allAutoChr/$chr
      for i in $window_bed.dxy; do
               cat $i | sed 's/$/ '$i'/' > $i.new
      done
      mv $window_bed.dxy.new ../dxy_allAutoChr/$chr
      rm $window_bed.dxy
done < ../pop_chr_info/$r_chr_bed

#### step 4 combine them all together
awk 'FNR == 3' $output/*.dxy.new > $output/$chr.100kWindow.waterbuck.dxy
sed 's/Pairwise neuclotide difference is:  //g' $output/$chr.100kWindow.waterbuck.dxy -i
cp $output/$chr.100kWindow.waterbuck.dxy ../
