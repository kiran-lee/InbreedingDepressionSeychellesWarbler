#!/bin/bash

#SBATCH --job-name=ROH
#SBATCH --output=ROH.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=4:00:00
#SBATCH --mail-user=kgllee1@sheffield.ac.uk
#SBATCH --mail-type=all
# #SBATCH -A molecolb
# #SBATCH -p molecolb

source /usr/local/extras/Genomics/.bashrc

SNPs=/fastdata/bop21kgl/RawData/allplates/FinalSTITCHImputationExtraSamples/mergedimputedchromosomeextrasamples.vcf.gz

#small ROH
plink \
  --allow-extra-chr \
  --vcf $SNPs \
  --geno 0.01 \
  --maf 0.01 \
  --homozyg-density 200 \
  --homozyg-gap 300 \
  --homozyg-het 2 \
  --homozyg-kb 375 \
  --homozyg-snp 50 \
  --homozyg-window-het 2 \
  --homozyg-window-missing 4 \
  --homozyg-window-snp 50 \
  --not-chr "ENA|OU383776|OU383776.1_RagTag" "ENA|OU383795|OU383795.1_RagTag" \
  --double-id \
  --out smallROHnew

#medium ROH
plink \
  --allow-extra-chr \
  --vcf $SNPs \
  --geno 0.01 \
  --maf 0.01 \
  --homozyg-density 200 \
  --homozyg-gap 300 \
  --homozyg-het 2 \
  --homozyg-kb 1360 \
  --homozyg-snp 50 \
  --homozyg-window-het 2 \
  --homozyg-window-missing 4 \
  --homozyg-window-snp 50 \
  --not-chr "ENA|OU383776|OU383776.1_RagTag" "ENA|OU383795|OU383795.1_RagTag" \
  --double-id \
  --out mediumROHnew

#large ROH
plink \
  --allow-extra-chr \
  --vcf $SNPs \
  --geno 0.01 \
  --maf 0.01 \
  --homozyg-density 200 \
  --homozyg-gap 300 \
  --homozyg-het 2 \
  --homozyg-kb 3750 \
  --homozyg-snp 50 \
  --homozyg-window-het 2 \
  --homozyg-window-missing 4 \
  --homozyg-window-snp 50 \
  --not-chr "ENA|OU383776|OU383776.1_RagTag" "ENA|OU383795|OU383795.1_RagTag" \
  --double-id \
  --out largeROHnew