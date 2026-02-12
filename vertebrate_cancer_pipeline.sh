#!/bin/bash
# Vertebrate Cancer Analysis Pipeline
# Comparative analysis of orthologous genes in birds and reptiles
# Detecting signatures of positive selection using HyPhy
# Originally developed for birds (Psittaciformes) and reptiles (Testudinata, Crocodylia)

set -e  # Stop on error

# USER-DEFINED PARAMETERS (you can pass them as arguments: ./script.sh SRC src Psittaciformes)
GENE_NAME="${1:-ATF3}"                 # e.g. SRC, ATF3
GENE_SYMBOL="${2:-atf3}"               # lowercase gene symbol e.g. atf3
GENE_ID="${3:-101939025}"             # only if needed (OPTIONAL) e.g. 101939025
ORTHOLOG_GROUP="${4:-testudinata,crocodylia}" # e.g. Psittaciformes or testudinata,crocodylia

# Configurable directories (change these if necessary)
BASE_DIR="${HOME}/vertebrate_cancer"  # Base directory (adjust to your preference)
WORKDIR="${BASE_DIR}/${GENE_NAME}_${ORTHOLOG_GROUP}"
RESULTS="${BASE_DIR}/Results/${GENE_NAME}_${ORTHOLOG_GROUP}"


# LOAD CONDA - Try multiple common paths if the first doesn't work
POSSIBLE_CONDA_PATHS=("${HOME}/miniconda3" "${HOME}/anaconda3" "/opt/miniconda3" "/opt/anaconda3")
CONDA_SOURCED=false

for CONDA_PATH in "${POSSIBLE_CONDA_PATHS[@]}"; do
    if [ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_PATH}/etc/profile.d/conda.sh"
        CONDA_SOURCED=true
        echo "Conda sourced from: ${CONDA_PATH}"
        break
    fi
done

if [ "$CONDA_SOURCED" = false ]; then
    echo "Error: Conda not found in common locations (${POSSIBLE_CONDA_PATHS[*]}). Install Miniconda/Anaconda and ensure it's in one of these paths, or add your path to POSSIBLE_CONDA_PATHS."
    exit 1
fi

# Create directories
mkdir -p "${WORKDIR}"
mkdir -p "${RESULTS}"
cd "${WORKDIR}"

echo "Starting pipeline for ${GENE_NAME} in ${ORTHOLOG_GROUP}..."

# 1. DOWNLOAD SEQUENCES
conda activate ncbi_datasets || { echo "Error: Could not activate ncbi_datasets environment."; exit 1; }

datasets download gene symbol "${GENE_SYMBOL}" \
  --include cds \
  --ortholog "${ORTHOLOG_GROUP}" \
  --filename "${GENE_NAME}.zip"

unzip "${GENE_NAME}.zip" -d data
tree data

conda deactivate

# Keep only CDS fasta files
find data/ncbi_dataset/data -type f ! -name "*.fna" -delete

# 2. HEADER CLEANING
cd data/ncbi_dataset/data
# Standardize headers
perl -pe 's/(>\w+_\d+\.\d)\:\d+\-\d+\s+.*organism\=(\w+)\s+(\w+)\$\s+.*/\1_\2\.\3/g' cds.fna > cds_clean.fna

# 3. ALIGNMENT WITH MUSCLE
if muscle -h 2>&1 | grep -q -- "-align"; then
  muscle -align cds_clean.fna -output "${GENE_NAME}_aligned.afa"
else
  muscle -in cds_clean.fna -out "${GENE_NAME}_aligned.afa"
fi

# 4. CLEAN CODONS WITH HYPHY
hyphy cln Universal "${GENE_NAME}_aligned.afa" No/Yes "${GENE_NAME}_aligned_clean.afa"

# 5. PHYLOGENETIC TREE (RAXML)
conda activate phylogenetics || { echo "Error: Could not activate phylogenetics environment."; exit 1; }

raxmlHPC \
  -s "${GENE_NAME}_aligned_clean.afa" \
  -n "${GENE_NAME}_raxml" \
  -m GTRCAT \
  -f a \
  -x 123 \
  -p 256 \
  -N autoMRE

conda deactivate

# Remove unnecessary RAxML files
find . -type f \( -name "RAxML_bootstrap*" -o -name "RAxML_bipartitions*" \) -delete

# 6. SELECTION ANALYSIS (HYPHY)
# Notes: Adjust LIBPATH if needed (e.g., to your Conda env)
LIBPATH="${CONDA_PATH}/lib/hyphy"

# aBSREL – branch selection
hyphy LIBPATH="${LIBPATH}" absrel \
  --alignment "${GENE_NAME}_aligned_clean.afa" \
  --tree "RAxML_bestTree.${GENE_NAME}_raxml" \
  --branches All \
  --output "${RESULTS}/${GENE_NAME}_aBSREL.json"

# MEME – site selection
hyphy LIBPATH="${LIBPATH}" meme \
  --alignment "${GENE_NAME}_aligned_clean.afa" \
  --tree "RAxML_bestTree.${GENE_NAME}_raxml" \
  --output "${RESULTS}/${GENE_NAME}_MEME.json"

# BUSTED – gene-wide selection
hyphy LIBPATH="${LIBPATH}" busted \
  --alignment "${GENE_NAME}_aligned_clean.afa" \
  --tree "RAxML_bestTree.${GENE_NAME}_raxml" \
  --branches All \
  --output "${RESULTS}/${GENE_NAME}_BUSTED.json"

echo "Pipeline completed successfully."
echo "Results saved in: ${RESULTS}"
echo "View results in HyPhy Vision: http://vision.hyphy.org"
