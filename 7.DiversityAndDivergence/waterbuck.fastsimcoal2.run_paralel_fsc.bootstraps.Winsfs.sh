#!/bin/bash

input=/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_NullModel_NonParametricBootstraps
obs_input=/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_NullModel_NonParametricBootstraps/get_InputFromWinsfs
batch=$1

for i in {1..100}
do

mkdir $input/waterbuck.TD.SamoleMatesis_boot_${batch}
mkdir $input/waterbuck.TD.SamoleMatesis_boot_${batch}/$i
cp $input/waterbuck.TD.SamoleMatesis.est $input/waterbuck.TD.SamoleMatesis_boot_${batch}/$i/$i.est
cp $input/waterbuck.TD.SamoleMatesis.tpl $input/waterbuck.TD.SamoleMatesis_boot_${batch}/$i/$i.tpl

export DIR=$input/waterbuck.TD.SamoleMatesis_boot_${batch}/$i
cd $DIR
cp $obs_input/$batch.2dsfs ${i}_jointDAFpop1_0.obs


/home/users/xi/software/fsc27_linux64/fsc2702 -t $i.tpl -n500000 -e $i.est -M -L100 -d -c10 -C 100 -r $i

done
