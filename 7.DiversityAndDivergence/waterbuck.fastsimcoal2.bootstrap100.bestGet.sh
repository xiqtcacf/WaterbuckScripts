input=/home/users/xi/Waterbuck_project/10.fastsimcoal27/fastsimcoal_NullModel_NonParametricBootstraps

batch=$1
export DIR=$input/waterbuck.TD.SamoleMatesis_boot_${batch}
cd $DIR

cp ../copy.sh ./
mkdir best
bash copy.sh $batch
cd best
for i in *.bestlhoods; do cat $i | sed 's/$/ '$i'/' > $i.new; done
awk 'FNR == 2' *.new > ../waterbuck.TD.SamoleMatesis.bootstrap.bestlhoods.50runs
cd ../
cp ../best_run.R ./
Rscript best_run.R
