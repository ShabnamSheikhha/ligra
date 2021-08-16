EXE_DIR=/home/shabsheikhha/Documents/BSC/ligra/apps
INPUT_DIR=/home/shabsheikhha/Documents/BSC/matrices/test
LOG_DIR=/home/shabsheikhha/Documents/BSC/ligra/apps/logs

MODE=$1
METHOD=$2
ALG=$3
DATASET=$4

EXE_FILE="$EXE_DIR/$ALG"
INPUT_FILE="$INPUT_DIR/adj-$DATASET.txt"
LOG_FILE="$LOG_DIR/$MODE-$METHOD-$DATASET.log"

command="$EXE_FILE $INPUT_FILE"
command_with_log="$EXE_FILE $INPUT_FILE > $LOG_FILE"

if [ $MODE = "PARALLEL" ] 
then
	export CILK=1
elif [ $MODE = "SEQUENTIAL" ] 
then 
	unset CILK
else
	echo "error"
	exit 1
fi

if [ $METHOD = "PARTITION" ]
then 
	export PARTITION=1
elif [ $METHOD = "NAIVE" ]
then 
	unset PARTITION
else
	echo "error"
	exit 1
fi

make clean > tmp.log
rm tmp.log
make -j $ALG
$command #> $LOG_FILE
