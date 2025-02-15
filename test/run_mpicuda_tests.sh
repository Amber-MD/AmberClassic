#!/bin/bash

if [ -z "$DO_PARALLEL" -o -z "$DO_CUDA" ]; then
   echo "Both DO_PARALLEL and DO_CUDA must be set to run mpicuda tests"
   exit 1
fi

source ../AmberClassic.sh 

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="$AMBERCLASSICHOME/logs"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.mpicuda.log"
difffile="${logprefix}.mpicuda.diff"

mkdir -p ${logdir}

echo "Running AmberClassic mpicuda tests on $date_string" | tee ${logfile}

(make --no-print-directory -k test.allcuda 2>&1) | tee -a ${logfile}

passed_count=`grep PASS ${logfile} | wc -l`
questionable_count=`grep "FAILURE:" ${logfile} | wc -l`
error_count=`grep "Program error" ${logfile} | grep -v "echo" | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
echo "${error_count} tests experienced errors" | tee -a ${logfile}

echo "Test log file saved as ${logfile}" | tee -a ${logfile}

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
else
   echo "No test diffs to save!" | tee -a ${logfile}
fi

if [ ${questionable_count} -gt 0 -o ${error_count} -gt 0 ]; then
    exit 1
fi
