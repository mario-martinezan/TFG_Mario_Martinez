#!/bin/bash

#SBATCH -c 32
#SBATCH --mem-per-cpu=3952M
#SBATCH -t 01:00:00
#SBATCH --gres=gpu:a100:1

./mdgpu_np setup.dat
