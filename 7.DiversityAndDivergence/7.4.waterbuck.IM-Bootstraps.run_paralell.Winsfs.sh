#!/bin/bash

input=/home/wlk579/0.waterbuck/fsc/0Using_boundry_nosingleton/bootstrap_IM
obs_input=/home/wlk579/0.waterbuck/fsc/0Using_boundry_nosingleton/bootstrap_IM/get_InputFromWinsfs
bootstrap=$1

for i in {1..20}
do

mkdir $input/waterbuck.IM.SamoleMatesis_boot_${bootstrap}
mkdir $input/waterbuck.IM.SamoleMatesis_boot_${bootstrap}/$i
cp $input/waterbuck.IM.SamoleMatesis.est $input/waterbuck.IM.SamoleMatesis_boot_${bootstrap}/$i/$i.est
cp $input/waterbuck.IM.SamoleMatesis.tpl $input/waterbuck.IM.SamoleMatesis_boot_${bootstrap}/$i/$i.tpl

export DIR=$input/waterbuck.IM.SamoleMatesis_boot_${bootstrap}/$i
cd $DIR
cp $obs_input/$bootstrap.2dsfs ${i}_jointDAFpop1_0.obs

/home/wlk579/Server_bos/apps/fsc27_linux64/fsc27093 -t $i.tpl \
-n500000 -e $i.est -M -L100 --nosingleton \
-d -c20 -C 100 -r $i

done
