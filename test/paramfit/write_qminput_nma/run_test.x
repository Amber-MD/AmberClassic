#!/bin/bash

echo "*** Create Gaussian Input Files (new prmtop) ***"

../../../bin/paramfit -i Job_Control.in -p mdcrd_calc/NMA.prmtop -c mdcrd_calc/mdcrd > prog_out.txt || { echo "  $0:  Program error"; exit 1; }

../../dacdif saved_output/Job.0.gjf.saved Job.0.gjf
../../dacdif saved_output/Job.1.gjf.saved Job.1.gjf
../../dacdif saved_output/Job.2.gjf.saved Job.2.gjf
../../dacdif saved_output/Job.3.gjf.saved Job.3.gjf
../../dacdif saved_output/Job.4.gjf.saved Job.4.gjf
/bin/rm prog_out.txt
