HgBenchmark_nested_NA

Developed by Yanxu Zhang
Last updated 13 Jan 2013

This directory contains a set of IDL scripts used to benchmark new and
updated versions of the GEOS-Chem mercury model at 0.5x0.666 North America grid. The primary script
used is MERCURY_BENCHMARK_NESTED_NA; this will call all other necessary routines.

Included in this directory are tracerinfo.dat and diaginfo.dat files.
However, you may want to replace these with the versions created by
your simulations.

The benchmarking scripts compare a "new" model run to a "reference"
simulation. If no reference file is specified, the new run will be
compared to the default file (included here) default.05x0666.NA.2009. 
This file uses the model version documented in Zhang et al. 2012 which is based on v9-01-02c.

To run the benchmarking procedures, enter the following at an IDL prompt:
   MERCURY_BENCHMARK_NESTED_NA, FILENAME=FILENAME, REFERENCE=REFERENCE, $
      PSFILENAME=PSFILENAME

Here, FILENAME is the name of the file containing output for the "new"
run, REFERENCE is the name of the file containing output for the
"reference" run, PSFILENAME is the name of the postscript in which the
output will be saved (default is mercury_benchmark_nested_na.ps).

If you have questions, please contact
Yanxu Zhang: zhangyanxv@gmail.com
