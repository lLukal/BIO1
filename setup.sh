#!/bin/bash

curl https://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta > ./input/MAP006-1_2D_pass.fasta

pip install -r ./src/requirements.txt

# Or
# conda install biopython

# To get started, use:
# bash run.sh --reference input/GCF_000005845.2_ASM584v2_genomic.fna --fragments input/MAP006-1_2D_pass.fasta --threads 4 --debug