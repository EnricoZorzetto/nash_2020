#!/bin/bash
echo 'Running nash analysis using all USHCN gauge data - REDO FAILED'
# n=1218
echo "Number of jobs [stations] = $n"
SECOND=$(sbatch -o stats.out -e stats.err --array=862 --parsable stats.q)
echo $SECOND
exit 0

