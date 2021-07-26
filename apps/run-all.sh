alg=$1
dataset=$2

echo "********************sequential naive********************"
source run.sh SEQUENTIAL NAIVE $alg $dataset
sleep 2.5

echo "********************sequential partitioned********************"
source run.sh SEQUENTIAL PARTITION $alg $dataset
sleep 2.5

echo "********************parallel naive********************" 
source run.sh PARALLEL NAIVE $alg $dataset
sleep 2.5 

echo "********************parallel partitioned********************" 
source run.sh PARALLEL PARTITION $alg $dataset

