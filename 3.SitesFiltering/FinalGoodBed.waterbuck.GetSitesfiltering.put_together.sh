batch=$1

ANGSD=/home/genis/software/angsd/angsd

all=/home/users/xi/Waterbuck_project/4.sites_filtering/6.all_together_twoREF/$batch/all.bed
sex_linked=/home/users/xi/Waterbuck_project/4.sites_filtering/6.all_together_twoREF/$batch/Ref.sexLinkedAndAbnormalScaff.txt
autosome=/home/users/xi/Waterbuck_project/4.sites_filtering/6.all_together_twoREF/$batch/nosexLinkedAndAbnormalScaff.bed
autosome_list=/home/users/xi/Waterbuck_project/4.sites_filtering/6.all_together_twoREF/$batch/nosexLinkedAndAbnormalScaff.list

#grep -vf $sex_linked $all | grep -v "mito" > $autosome
grep -vf $sex_linked $all > $autosome
#cut -f1 $autosome > $autosome_list
cut -f1 $autosome > $autosome_list

BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools
outdir=/home/users/xi/Waterbuck_project/4.sites_filtering/6.all_together_twoREF/$batch


rep=$outdir/Ref.Nonrepeat.bed
map=$outdir/Ref.mappability_m1_k150_e2.bed
dep=$outdir/Depth.all_keep.bed
het=$outdir/Excess_heter.good.bed


$BEDTOOLS intersect -a $autosome -b $rep > $outdir/nosexLinkedAndAbnormalScaff_rep.bed 
$BEDTOOLS intersect -a $outdir/nosexLinkedAndAbnormalScaff_rep.bed -b $het > $outdir/nosexLinkedAndAbnormalScaff_rep_het.bed
$BEDTOOLS intersect -a $outdir/nosexLinkedAndAbnormalScaff_rep_het.bed -b $dep > $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep.bed
$BEDTOOLS intersect -a $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep.bed -b $map > $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep_map.bed

awk '{print $1"\t"$2+1"\t"$3}' $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep_map.bed > $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions

$ANGSD sites index $outdir/nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions
