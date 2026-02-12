# Vertebrate Cancer Analysis Pipeline

## Overview
This repository contains a **jointly developed Linux-based bioinformatics pipeline**
for the **comparative analysis of orthologous genes across vertebrates**, with a
focus on **birds (Psittaciformes)** and **reptiles (Testudinata, Crocodylia)**.

The pipeline detects **signatures of positive selection** potentially associated
with cancer-related traits using codon-aware evolutionary models.


## Pipeline workflow
The pipeline performs the following steps:

1. Download coding sequences (CDS) of orthologous genes from **NCBI datasets**
2. Perform multiple sequence alignment using **MUSCLE**
3. Clean FASTA headers and remove stop codons using **HyPhy clean**
4. Infer maximum-likelihood phylogenetic trees using **RAxML**
5. Detect signatures of selection using **HyPhy**:
   - aBSREL (branch-level selection)
   - MEME (site-level selection)
   - BUSTED (gene-wide episodic selection)

All steps are executed in a **Linux terminal environment**.


## Requirements
- Linux OS
- Conda
- NCBI datasets / Entrez Direct
- MUSCLE
- RAxML
- HyPhy
- iTOL (optional, tree visualization)


## Usage

Edit the following variables at the top of the script:

- `GENE_NAME, GENE_SYMBOL` – Gene details (e.g. SRC) 
- `ORTHOLOG_GROUP` – Taxonomic group (birds or reptiles) eg. Psittaciformes, Testudinata, or Crocodylia

Run the pipeline:

```bash
bash vertebrate_cancer_pipeline.sh
```

## Outputs
- Codon-aware multiple sequence alignments
- Maximum-likelihood phylogenetic trees
- HyPhy JSON output files:
- aBSREL (branches under episodic positive selection)
- MEME (sites under positive selection)
- BUSTED (gene-wide episodic positive selection across specified lineages)

Results can be visualized using HyPhy Vision:
http://vision.hyphy.org

## Notes
- The pipeline is designed for Linux terminal usage.
- Sequence headers are cleaned to ensure compatibility with downstream analyses.
- The workflow can be adapted to analyze other genes or vertebrate taxa.
- Output files are suitable for downstream visualization and statistical analysis.