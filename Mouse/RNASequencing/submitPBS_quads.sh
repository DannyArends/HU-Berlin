#!/bin/bash
#PBS -N extractExonExpression
#PBS -l nodes=1:ppn=6
#PBS -q quads
#PBS -l walltime=20:00:00
#PBS -lpmem=16384M

cd /data/p256802/
Rscript extractExonExpression.R
