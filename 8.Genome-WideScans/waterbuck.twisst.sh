# YOU NEED TO INCLUDE REFERENCE & SITES PATHS
# You need a bam list file with path to each bam file
# TWISST GROUP NAMES ARE HARD-CODED
# groupFile needs to be made separately knowing number of samples and numbered as the order of the bam list file. Second column is population name.
## for ((i=1;i<=46;i++)); do sed -n ${i}p bams_popinfo_twisstv1.txt | cut -f1 | awk -v num="$i" '{print "V"num"\t"$1}' >> groups.txt ; done
# Before starting, set-up a temp folder (or add it to script eventually)

## To use e.g. bash run_twisst.sh 100000 CHR_01 > something_CHR_01.out 2>&1

REF="/brendaData/xi/African1Kg/Waterbuck/Reference/Goat.fasta"
SITES="/home/users/xi/Waterbuck_project/5.genotype_likelihood/Goat.29Chr.sites_filter/Goat.29Chr.nosexLinkedAndAbnormalScaff_rep_het_dep_map.regions"
window=$1
contig=$2
chr=$3
win=$(echo $window/1000 | bc)
size=$(grep -w ${contig} ${REF}.fai | cut -f2)
start=${window}

# This takes the window size of interest and divides the scaffold into those windows
while [ $start -le ${size} ]
do
 stop=$(($start + ${window}))
 echo "$contig:$start-$stop" >> regions_${win}k_${contig}.txt
 start=$stop
done

# This feeds the above regions file into angsd to obtain an IBSmatrix followed by an Rscript which creates a NJ tree
## THIS PART NEEDS TO BE REPLACED WHEN USING PHASED DATA
cat regions_${win}k_${contig}.txt | while read line
do
 /home/users/xi/software/angsd/angsd -GL 2 -out temp/${line} -nThreads 10 -bam bams_5pops.txt -minmapQ 30 -minQ 20 -r $line -sites ${SITES} -doIBS -1 -makeMatrix 1 -doCounts 1
 cd temp
 Rscript ../plotNJtree_regions.R ${line}.ibsMat ${contig}_${win}k_njTree.txt
 rm ${line}.*
 cd ../
done

# CLEAN UP FILES RELATED TO WINDOWS WHERE NO NJ TREE WAS PRODUCED; FOR PLOTTING PURPOSES.

#paste regions_${win}k_${contig}.txt <(grep "filtering" something_${contig}.out | cut -d ":" -f2 | sed 's/ //g') > regions_withNumSites_${win}k_${contig}.txt
#grep -B7 "Error" something_${contig}.out | grep "temp" | cut -d/ -f2 | cut -d. -f1 |cut -d: -f2 > ${contig}_errors.txt
#grep -w -v -f ${contig}_errors.txt regions_withNumSites_${win}k_${contig}.txt | sed -e 's/:/\t/g' -e 's/-/\t/g' - > clean_regions_withNumSites_${win}k_${contig}.txt
#rm regions_withNumSites_${win}k_${contig}.txt
#echo -e "chrom\tstart\tend\tsites" | cat - clean_regions_withNumSites_${win}k_${contig}.txt > tmp_${contig}
#mv tmp_${contig} clean_regions_withNumSites_${win}k_${contig}.txt

#xi rewrite this part:
paste regions_${win}k_${contig}.txt <(grep "filtering" waterbuck.${win}k.${contig}.${chr}.out | cut -d ":" -f2 | sed 's/ //g') > regions_withNumSites_${win}k_${contig}.txt
grep -B7 "Error" waterbuck.${win}k.${contig}.${chr}.out | grep "temp" | cut -d/ -f2 | sed 's/.ibsMat"//g' > ${win}k.${contig}_errors.txt
grep -w -v -f ${win}k.${contig}_errors.txt regions_withNumSites_${win}k_${contig}.txt | sed -e 's/:/\t/g' -e 's/-/\t/g' - > clean_regions_withNumSites_${win}k_${contig}.txt
rm regions_withNumSites_${win}k_${contig}.txt
echo -e "chrom\tstart\tend\tsites" | cat - clean_regions_withNumSites_${win}k_${contig}.txt > tmp_${contig}.${win}k
mv tmp_${contig}.${win}k clean_regions_withNumSites_${win}k_${contig}.txt

# RUN TWISST: inputs are trees, groupsFile, and group names [-g]. THESE ARE HARD-CODED.
python3 ~/software/twisst/twisst.py -t temp/${contig}_${win}k_njTree.txt_NJtrees.txt -w ${contig}_${win}k.weights.csv -g Kafue -g Samole -g CentralDef -g NorthCom -g SouthCom --groupsFile groups.txt --method complete --outputTopos topos_${contig}_${win}k.txt
