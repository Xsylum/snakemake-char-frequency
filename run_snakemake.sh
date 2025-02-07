#!/bin/bash
#SBATCH --job-name=text_count_pipeline
#SBATCH --partition=standard
#SBATCH --account=aafc_aac
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000M
#SBATCH --comment="registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu20.04"

# PRECONDITION: Ensure you have changed your initial conda environment to the one provided in this repository

snakemake --profile profile