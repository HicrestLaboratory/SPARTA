#!/bin/bash

for Matrix in data/minitest/*.el
do
    sbatch batch/csrMult $Matrix 16
    sbatch batch/BellFix $Matrix 16
    sbatch batch/BellMult $Matrix 16
    sbatch batch/VbrFix $Matrix 16
    sbatch batch/VbrNoFix $Matrix 16
done

for Matrix in data/minitest/*
do
    sbatch batch/csrMult $Matrix 32
    sbatch batch/BellFix $Matrix 32
    sbatch batch/BellMult $Matrix 32
    sbatch batch/VbrFix $Matrix 32
    sbatch batch/VbrNoFix $Matrix 32
done

for Matrix in data/minitest/*
do
    sbatch batch/csrMult $Matrix 64
    sbatch batch/BellFix $Matrix 64
    sbatch batch/BellMult $Matrix 64
    sbatch batch/VbrFix $Matrix 64
    sbatch batch/VbrNoFix $Matrix 64
done
