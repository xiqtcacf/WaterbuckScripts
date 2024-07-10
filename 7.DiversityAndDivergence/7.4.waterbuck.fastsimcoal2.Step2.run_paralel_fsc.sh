#!/bin/bash

batch=$1  ###model name: NULL, IM, MStop
input="/home/wlk579/0.waterbuck/fsc/0Using_boundry_nosingleton/$batch"

cd ${input}

for i in {1..100}; do

cp ../$batch.SamoleMatesis.est $batch.SamoleMatesis_${i}.est
cp ../$batch.SamoleMatesis.tpl $batch.SamoleMatesis_${i}.tpl
cp ../_jointDAFpop1_0.obs $batch.SamoleMatesis_${i}_jointDAFpop1_0.obs

/home/wlk579/Server_bos/apps/fsc27_linux64/fsc27093 -t $batch.SamoleMatesis_${i}.tpl \
-n500000 -e $batch.SamoleMatesis_${i}.est -M -L100 --nosingleton \
-d -c10 -C 100 -r $i

done
