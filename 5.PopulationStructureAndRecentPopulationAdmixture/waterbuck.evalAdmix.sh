#!/bin/bash

batch=$1
export DIR=/home/users/xi/Waterbuck_project/6.population_structure/NGSadmix/best_Xi/
cd $DIR

ngsadmix=/home/users/xi/software/ngsadmix/NGSadmix
evalAdmix=/home/users/xi/software/evalAdmix/evalAdmix
beagle=/home/users/xi/Waterbuck_project/6.population_structure/waterbuck.$batch.all.beagle.gz

###### evaluate the results of an admixture analysis

$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K10.3.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K10.3.qopt -o evaladmixOut.$batch.K10.3.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K2.39.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K2.39.qopt -o evaladmixOut.$batch.K2.39.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K3.37.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K3.37.qopt -o evaladmixOut.$batch.K3.37.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K4.45.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K4.45.qopt -o evaladmixOut.$batch.K4.45.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K5.10.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K5.10.qopt -o evaladmixOut.$batch.K5.10.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K6.18.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K6.18.qopt -o evaladmixOut.$batch.K6.18.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K7.58.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K7.58.qopt -o evaladmixOut.$batch.K7.58.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K8.14.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K8.14.qopt -o evaladmixOut.$batch.K8.14.corres -P 20
$evalAdmix -beagle $beagle -fname waterbuck.DefassaWaterbuckRef.K9.42.fopt.gz -qname waterbuck.DefassaWaterbuckRef.K9.42.qopt -o evaladmixOut.$batch.K9.42.corres -P 20
