#!/bin/bash

# Defines some color codes to use
RED='\033[0;31m'
BLUE='\033[0;34m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
NOCOLOR='\033[0m'

# Parse the command line. We need 2 arguments: the language to benchmark
# and the number of threads to use

show_help() 
{
    echo "USAGE: run_benchmark [--lang={C,CHAPEL}] [--threads=NUMTHREADS]"
}

# default values
NUMTHREADS=1
LANGUAGE="CHAPEL"

# parse arguments
for i in "${@}"; do
    case "${i}" in
        # help
        -h|--help)
        show_help
        exit 0
        ;;
        --lang=*)
        LANGUAGE=${i#*=}
        ;;
        --threads=*)
        NUMTHREADS=${i#*=}
        ;;
        # bad argument
        *)
        die "Unknown option '${i}'"
        ;;
    esac
done

# Store language in all uppercase
LANGUAGE=${LANGUAGE^^}

# Make sure language is valid
if [[ "${LANGUAGE}" != "CHAPEL" && "${LANGUAGE}" != "C" ]];
then
    printf "${RED} ERROR: Language specified with --lang must be C or Chapel (case insensitive)\n${NOCOLOR}"
    exit
fi

# Common directories to refer to
HOMEDIR=/home/tbrolin
THISDIR=$PWD
TENSORDATA=${HOMEDIR}/tensor_data/CHAPEL
OUTDIR=${THISDIR}/output_data_${LANGUAGE}

# Data sets
declare -a dataSets=("YELP" "NELL-2")

######################################################################################
#
# 1.) Create output directory
#
######################################################################################
echo ""
printf "${YELLOW}######################################################################################\n${NOCOLOR}"
totalStart=`date +%s`
currDate=`date`
printf "${RED}+ START TIME: ${currDate}\n\n${NOCOLOR}"

printf "${CYAN}+ Creating output directory: $OUTDIR\n${NOCOLOR}"

mkdir -p ${OUTDIR}
for dataSet in "${dataSets[@]}"
do
    mkdir -p ${OUTDIR}/${dataSet}/${NUMTHREADS}THD
done

#######################################################################################
##
## 2.) Run trials
##
#######################################################################################
echo ""
EXEPATH=${THISDIR}/execs
TIME_START=`date +%s`
currDate=`date`
printf "${YELLOW}######################################################################################\n${NOCOLOR}"
printf "${CYAN}+ Running tests for ${LANGUAGE} on ${NUMTHREADS} threads\n${NOCOLOR}"
printf "${RED} + START TIME: ${currDate}\n${NOCOLOR}"

for dataSet in "${dataSets[@]}"
do
    DATAPATH="$TENSORDATA/${dataSet}/*.bin"
    printf "${GREEN}\t- Running on ${dataSet} data set\n${NOCOLOR}"
    OUTFILEDIR=${OUTDIR}/${dataSet}/${NUMTHREADS}THD
    OPTS=""
    if [ ${LANGUAGE} == "C" ];
    then
        OPTS="cpd ${DATAPATH} --rank=35 --seed=123456 --tol=1e-20 -i 20 --threads=${NUMTHREADS} --nowrite -vv"
    else
        OPTS="${DATAPATH} --rank=35 --seed=123456 --tol=1e-20 --iters=20 --threads=${NUMTHREADS}"
    fi
    for trial in 1 2 3 4 5
    do
        currDate=`date`
        OUTFILENAME="${LANGUAGE}_${dataSet}_${NUMTHREADS}THD_Trial${trial}.txt"
        printf "${BLUE}\t\t- [$currDate] Data Set=${dataSet}; Threads=$NUMTHREADS; Trial=$trial\n${NOCOLOR}"
        printf "${BLUE}\t\t\tOMP_NUM_THREADS=${NUMTHREADS} ${EXEPATH}/./splatt_${LANGUAGE} ${OPTS}\n"
        OMP_NUM_THREADS=${NUMTHREADS} ${EXEPATH}/./splatt_${LANGUAGE} ${OPTS} > ${OUTFILEDIR}/${OUTFILENAME}
    done
    printf "${GREEN}\t- Done\n${NOCOLOR}"

TIME_END=`date +%s`
TIME_TOTAL=$((TIME_END-TIME_START))
currDate=`date`
printf "${RED}END TIME:   ${currDate}\n${NOCOLOR}"
printf "${RED}TOTAL TIME: ${TIME_TOTAL} seconds\n\n${NOCOLOR}"
done


##########################################################################
totalEnd=`date +%s`
totalTime=$((totalEnd-totalStart))
currDate=`date`
echo ""
printf "${YELLOW}######################################################################################\n${NOCOLOR}"
printf "${RED}END TIME:   ${currDate}\n${NOCOLOR}"
printf "${RED}TOTAL TIME: ${totalTime} seconds\n${NOCOLOR}"
echo ""
