# README

## This repository is a compilation of scripts from Lareau Lab, Ott Lab, and Fletcher Lab contributions. It has been specifically modified for Cas12 compatibility (while retaining Cas13 functionality)
* https://github.com/lareaulab/cas13a_guide_design
* https://github.com/duopeng/Cas13a_guides_Influenza
* (will enter Ott lab link here)

## Data needed
reference genome sequence (.fa)
sequenced viral genomes representing genetic diversity (.fa)
Cas13a repeat sequence

## Need to install:
bowtie
RNAfold
R, with libraries: optparse, here, Biostrings, parallel, doParallel
Workflow (using scripts)

## Documentation
1. Pull genomes from NCBI, etc. Extract metadata and fasta.
   Script: extract_metadata_and_fasta.py
   This script needs to be modified for query and filtering on lines 35, 85-100, 109, and 112.
2. Split fasta file by subtype
   Script: split_fasta_by_subtype.py
   Inputs: metadata file, fasta sequences
   Bash command: python coreScripts/split_fasta_bySubtype.py --metadata Amazonensis_data/Leishmania_amazonensis_metadata.tsv --fasta Amazonensis_data/Leishmania_amazonensis_sequences.fasta
3. Break each genome into 20nt windows
   Script: break_genome.R
   Inputs: reference genome .fa file, enzyme name, output directory
   Bash command: python coreScripts/break_genomes.py --enzyme Cas12 --scriptdir coreScripts/ --input subtype_Leishmania_amazonensis.fa --out tiled_genomes/tiled_subtype_amazonensis --num_workers 32
4. Filter for inclusivity of guides threshold (eg. 50% of strains hit)
   Script: reformatFilterGuideHits.R
   Input: threshold number
   Bash command: Rscript coreScripts/reformatFilterGuideHits.R leishmania_all .95
5. Calculate crRNA folding structures
   Script: score_RNAgold_crRNAs.R
   Input: enzyme name, output directory, RNAfold, cas repeat seqeunce (Cas12 sequences show in example below)
   Bash command: Rscript coreScripts/score_RNAfold_crRNAs.R --enzyme Cas12 --rnafold RNAfold --out . --cas_repeat UAAUUUCUACUAAGUGUAGAU
6. Filter guides for Hairpin formation and fold score
   Script: filterGuidesFoldScore.R
   Bash command: Rscript coreScripts/filterGuidesFoldScore.R leishmania_all

At this point you will have a list of guide candidates that you can now filter check for cross-reactivity. There are many ways you cna do this (bowtie, Kraken, etc.), but I would recommmend putting the guides through NCBI BLAST and filtering for only those that have less than 1/2 mismatches and seqeunce coverage. Remember if you are using Cas12, the PAM sequences needs to be included on the reverse complement region and needs to match
