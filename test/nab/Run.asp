#!/bin/sh

. ../program_error.sh

if [ "$DO_PARALLEL" != "" ];  then
    NAB="mpinab"
else
    NAB="nab"
fi

echo "Running test to compute GB Newton-Raphson and normal modes:"
echo ""
../../bin/$NAB -o asp asp.nab || error
($DO_PARALLEL ./asp > asp.out.tmp) || error
cat asp.out.tmp | tail -67 > asp.out

echo "Note: Very small differences between asp.out and the saved file are"
echo "common and do not necessarily indicate a problem."
../dacdif asp.out.check asp.out

rm -f asp asp.c asp.out.tmp
exit 0
