#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -J cori_bigmem
#SBATCH -t 00:20:00
#SBATCH --mem=700G
#SBATCH -M escori
#SBATCH --array=1-5
#SBATCH -p bigmem
set -e;

bash array.bash
