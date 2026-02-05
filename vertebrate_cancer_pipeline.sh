#!/bin/bash
# Vertebrate Cancer Analysis Pipeline
# Comparative analysis of orthologous genes in birds and reptiles
# Detecting signatures of positive selection using HyPhy
# Originally developed for birds (Psittaciformes) and reptiles (Testudinata, Crocodylia) 

# USER-DEFINED PARAMETERS 

# Gene information
GENE_NAME="SRC"                 # e.g. SRC, ATF3
GENE_SYMBOL="src"               # lowercase gene symbol
GENE_ID="101939025"             # only used if needed

# Taxa options:
# Birds: Psittaciformes
# Reptiles: testudinata,crocodylia
ORTHOLOG_GROUP="Psittaciformes"   # change to testudinata,crocodylia for reptiles

# Working directories
WORKDIR="$HOME/vertebrate_cancer/${GENE_NAME}_${ORTHOLOG_GROUP}"
RESULTS="$HOME/Results/${GENE_NAME}_${ORTHOLOG_GROUP}"

# LOAD CONDA
source /home/manager/anaconda3/etc/profile.d/conda.sh

mkdir -p ${WORKDIR}
mkdir -p ${RESULTS}
cd ${WORKDIR}

# 1. DOWNLOAD SEQUENCES
conda activate ncbi_datasets

# Option A: download by gene symbol (used for birds)
datasets download gene symbol ${GENE_SYMBOL} \
  --include cds \
  --ortholog ${ORTHOLOG_GROUP} \
  --filename ${GENE_NAME}.zip

# Option B: download by gene ID (used if needed)
# Uncomment if needed
# datasets download gene gene-id ${GENE_ID} \
#   --include cds \
#   --ortholog ${ORTHOLOG_GROUP} \
#   --filename ${GENE_NAME}.zip

unzip ${GENE_NAME}.zip -d data
tree data

conda deactivate

# Keep only CDS fasta files
find data/ncbi_dataset/data -type f ! -name "*.fna" -delete

# 2. HEADER CLEANING
cd data/ncbi_dataset/data
# Standardize headers
perl -pe 's/(>\w+_\d+\.\d)\:\d+\-\d+\s+.*organism\=(\w+)\s+(\w+)\]\s+.*/\1_\2\.\3/g' cds.fna > cds_clean.fna

# 3. ALIGNMENT WITH MUSCLE
muscle -in cds_clean.fna -out ${GENE_NAME}_aligned.afa

# 4. CLEAN CODONS WITH HYPHY
hyphy cln Universal ${GENE_NAME}_aligned.afa No/Yes ${GENE_NAME}_aligned_clean.afa

# 5. PHYLOGENETIC TREE (RAXML)

conda activate phylogenetics

raxmlHPC \
  -s ${GENE_NAME}_aligned_clean.afa \
  -n ${GENE_NAME}_raxml \
  -m GTRCAT \
  -f a \
  -x 123 \
  -p 256 \
  -N autoMRE

conda deactivate

# Remove unnecessary RAxML files
find . -type f \( -name "RAxML_bootstrap*" -o -name "RAxML_bipartitions*" \) -delete
# 6. SELECTION ANALYSIS (HYPHY)

# Notes:
# - LIBPATH specifies HyPhy library directory
# - Alignment file: codon alignment (from MUSCLE)
# - Tree file: best phylogenetic tree (from RAxML)
# - Output: results saved in JSON format

# aBSREL – branch selection
hyphy LIBPATH=/home/manager/anaconda3/lib/hyphy absrel \
  --alignment ${GENE_NAME}_aligned_clean.afa \
  --tree RAxML_bestTree.${GENE_NAME}_raxml \
  --branches All \
  --output ${RESULTS}/${GENE_NAME}_aBSREL.json

# MEME – site selection
hyphy LIBPATH=/home/manager/anaconda3/lib/hyphy meme \
  --alignment ${GENE_NAME}_aligned_clean.afa \
  --tree RAxML_bestTree.${GENE_NAME}_raxml \
  --output ${RESULTS}/${GENE_NAME}_MEME.json

# BUSTED – gene-wide selection (used in reptiles)
hyphy LIBPATH=/home/manager/anaconda3/lib/hyphy busted \
  --alignment ${GENE_NAME}_aligned_clean.afa \
  --tree RAxML_bestTree.${GENE_NAME}_raxml \
  --branches All \
  --output ${RESULTS}/${GENE_NAME}_BUSTED.json

# View results in HyPhy Vision (http://vision.hyphy.org)

echo "Pipeline completed successfully."
echo "Results saved in: ${RESULTS}"
