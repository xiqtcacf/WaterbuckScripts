#!/bin/bash

batch=$1
export DIR=/home/users/xi/Waterbuck_project/6.population_structure/PCAngsd/
cd $DIR

pcangsd=/home/users/xi/software/pcangsd/pcangsd.py
beagle=/home/users/xi/Waterbuck_project/6.population_structure/waterbuck.$batch.all.beagle.gz

for i in {1..10}
do
python3 $pcangsd -beagle $beagle -e $i -o waterbuck.$batch.Pcangsd.PC${i} -minMaf 0.05 -threads 20 
done
