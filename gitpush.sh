#!/bin/sh

#SBATCH --time=01-00:00:00

git add --all
git commit -m "reduce size"
