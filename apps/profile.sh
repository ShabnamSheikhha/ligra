algorithm=$1
DEBUG=false
MODE=SEQUENTIAL
METHOD=PARTITION

unset CILK
#export CILK=1
export PARTITION=1
#unset PARTITION
make clean
make -j $algorithm

MATRIX_DIR=/home/shabsheikhha/Documents/BSC/matrices/simple
PERF_DIR=/home/shabsheikhha/Documents/BSC/ligra/apps/perf
ALGORITHM_DIR=/home/shabsheikhha/Documents/BSC/ligra/apps

#datasets=(amazon-2008 cnr-2000 eu-2005 hollywood-2009 in-2004)
datasets=(simple-1024-4096)

general=cycles,instructions,cache-references,cache-misses,branches,branch-misses,bus-cycles,context-switches
L1=L1-dcache-stores,L1-dcache-loads,L1-dcache-load-misses
LLC=LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses

for dataset in "${datasets[@]}"
do 
	echo "perf for $dataset"
	
	perf_file="$PERF_DIR/$algorithm-$MODE-$METHOD-$dataset.data"
	matrix_file="$MATRIX_DIR/adj-$dataset.txt"
	algorithm_exe="$ALGORITHM_DIR/$algorithm"

	record="sudo perf record -e $general,$L1,$LLC -o $perf_file -a  $algorithm_exe $matrix_file"
	report="sudo perf report -i $perf_file"
	
	echo $record
	$record
	echo "run '$report' to see results"
	echo "*************************"
	echo "*************************"
	sleep 5
done

#echo $command
