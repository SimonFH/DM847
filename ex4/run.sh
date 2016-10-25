#!/bin/bash
if [ $# -eq 0 ]; then
    rscript SEQMOTIF.R ex4_upstream.fas.txt 19 && rscript selo.R candidates.csv -t DNA -s IC -f CSV -bg bgdist.csv
elif [ $1 = "clean" ]; then
    find . -name "out.png"\
    -or -name "consensus.txt"\
    -or -name "candidates.csv"\
    -or -name "bgdist.csv"\
    -or -name "PCM.csv"\
    -or -name "PWM.csv"|xargs -I % rm %
fi
