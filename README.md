#hi

python coreScripts/split_fasta_bySubtype.py --metadata Amazonensis_data/Leishmania_amazonensis_metadata.tsv --fasta Amazonensis_data/Leishmania_amazonensis_sequences.fasta

python coreScripts/break_genomes.py --enzyme Cas12 --scriptdir coreScripts/ --input subtype_Leishmania_amazonensis.fa --out tiled_genomes/tiled_subtype_amazonensis --num_workers 32

python coreScripts/count_combine_guides.py --metadata Amazonensis_data/Leishmania_amazonensis_metadata.tsv --pathPrefix tiled_genomes/tiled_subtype_

Rscript coreScripts/reformatFilterGuideHits.R leishmania_all .95

Rscript coreScripts/score_RNAfold_crRNAs.R --enzyme Cas12 --rnafold RNAfold --out . --cas_repeat UAAUUUCUACUAAGUGUAGAU

Rscript coreScripts/filterGuidesFoldScore.R leishmania_all
