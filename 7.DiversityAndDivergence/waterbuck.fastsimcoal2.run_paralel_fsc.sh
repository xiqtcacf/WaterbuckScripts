
#!/bin/bash

num_procs=7
num_jobs="\j"


input="/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_exe"


cd ${input}

for i in {1..100}; do

   while (( ${num_jobs@P} >= num_procs )); do
      #echo "${num_jobs@P} "
    wait -n
   done

cp waterbuck.TD.SamoleMatesis.est waterbuck.TD.SamoleMatesis_${i}.est
cp waterbuck.TD.SamoleMatesis.tpl waterbuck.TD.SamoleMatesis_${i}.tpl
cp waterbuck.TD.SamoleMatesis_jointDAFpop1_0.obs waterbuck.TD.SamoleMatesis_${i}_jointDAFpop1_0.obs

nice -n10 /home/users/xi/software/fsc27_linux64/fsc2702 -t waterbuck.TD.SamoleMatesis_${i}.tpl \
-n500000 -e waterbuck.TD.SamoleMatesis_${i}.est -M -L100 \
-d -c10 -C 100 -r ${i} 2> ${input}/waterbuck.TD.SamoleMatesis_${i}.log &

done
