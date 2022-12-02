#!/bin/bash

#file=/home/users/xi/Waterbuck_project/6.population_structure/waterbuck.DefassaWaterbuckRef.all.beagle.gz
nfile=waterbuck.DefassaWaterbuckRef
num=100  # Maximum number of iterations


K=$1 #number of populations. 
out=/home/users/xi/Waterbuck_project/6.population_structure/NGSadmix/$K #output  directory (no /)
out2=/brendaData/xi/African1Kg/Waterbuck/6.population_structure/NGSadmix

for f in {1..100}
do
    grep "like=" $out/$nfile.$K.$f.log | cut -f2 -d " " | cut -f2 -d "=" >> $out/$nfile.$K.likes
    CONV=`Rscript -e "r<-read.table('$out/$nfile.$K.likes');r<-r[order(-r[,1]),];cat(sum(r[1]-r<3),'\n')"` #Check for convergence. 3,5...?
    if [ $CONV -gt 2 ]  #-gt 2 = greater than 2
    then
	cp $out/$nfile.$K.$f.qopt $out2/$nfile.$K.$f.qopt_conv
	cp $out/$nfile.$K.$f.fopt.gz $out2/$nfile.$K.$f.fopt_conv.gz
	cp $out/$nfile.$K.$f.log $out2/$nfile.$K.$f.log_conv
	break
    fi
done

