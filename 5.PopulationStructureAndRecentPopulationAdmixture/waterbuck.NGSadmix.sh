#!/bin/bash

batch1=$1
batch2=$2
export DIR=/home/users/xi/Waterbuck_project/6.population_structure/NGSadmix
cd $DIR

ngsadmix=/home/users/xi/software/ngsadmix/NGSadmix
beagle=/home/users/xi/Waterbuck_project/6.population_structure/waterbuck.$batch1.all.beagle.gz

for i in {2..100}
do
$ngsadmix -likes $beagle -K $batch2 -o K${batch2}/waterbuck.$batch1.K${batch2}.$i -minMaf 0.05 -P 5
done
