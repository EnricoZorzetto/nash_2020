#!/bin/bash
mkdir -p ../Reports
echo 'Running nash analysis using all USHCN gauge data'
n=1218
echo "Number of jobs [stations] = $n"
SECOND=$(sbatch -o stats.out -e stats.err --array 1-$n --parsable stats.q)
echo $SECOND
exit 0

