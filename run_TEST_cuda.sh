srun --job-name=TEST_cuda --ntasks=1 --partition=short --gres=gpu:1 ./programs/cuda/TEST_cuda -b 3 -v 2 -a 2 -B 3
