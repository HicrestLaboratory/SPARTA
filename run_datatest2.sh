#!/bin/bash

path=../SPARTA_datasets/suitsparse_collection_3_no3-5

echo "sbatch batch/VBR $path/ckt11752_dc_1.mtx 128 512"
sbatch batch/VBR_rect $path/ckt11752_dc_1.mtx 128 512 0.2

echo "sbatch batch/VBR $path/email-Enron.mtx 128 512"
sbatch batch/VBR_rect $path/email-Enron.mtx 128 512 0.05

echo "sbatch batch/VBR $path/FEM_3D_thermal1.mtx 128 512"
sbatch batch/VBR_rect $path/FEM_3D_thermal1.mtx 128 512 0.01

echo "sbatch batch/VBR $path/g7jac100sc.mtx 128 512"
sbatch batch/VBR_rect $path/g7jac100sc.mtx 128 512 0.01

echo "sbatch batch/VBR $path/k3plates.mtx 128 512"
sbatch batch/VBR_rect $path/k3plates.mtx 128 512 0.09

echo "sbatch batch/VBR $path/lhr17.mtx 128 512"
sbatch batch/VBR_rect $path/lhr17.mtx 128 512 0.08

echo "sbatch batch/VBR $path/TEM27623.mtx 128 512"
sbatch batch/VBR_rect $path/TEM27623.mtx 128 512 0.05

for Matrix in ckt11752_dc_1.mtx email-Enron.mtx FEM_3D_thermal1.mtx g7jac100sc.mtx k3plates.mtx lhr17.mtx TEM27623.mtx
do
    echo "sbatch batch/CSR $path/$Matrix"
    sbatch batch/CSR $path/$Matrix
done
