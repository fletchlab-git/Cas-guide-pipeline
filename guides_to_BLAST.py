#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:39:15 2025

@author: zainamoussa
"""

from Bio.Seq import Seq
import csv


def process_sequences(input_txt, output_fasta):
    with open(input_txt, 'r') as infile, open(output_fasta, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')  # Assuming tab-delimited file
        for row in reader:
            accession = row['accession'].strip()
            target = row['target'].strip()
            pam = row['PAM'].strip()
            strand = row['strand'].strip()
            start = row['start'].strip()
            
            target_seq = Seq(target)
            pam_seq = Seq(pam)
            
            if strand == '+':
                result_seq = target_seq + pam_seq.reverse_complement()
            elif strand == '-':
                result_seq = pam_seq + target_seq.reverse_complement()
            else:
                continue  # Skip rows with invalid strand information
            print(result_seq)
            outfile.write(f">{accession}_{start}\n{result_seq}\n")




# Example usage
process_sequences('/Users/zainamoussa/Desktop/Leishmania_viannia_pipeline/top_guide_crossreactivity/Cas12_cross_reactivity_guide_results_summary copy.txt', 'BLAST_viannia_guides.fa')

#%%
def split_fasta(input_fasta, output_prefix, num_files=5):
    with open(input_fasta, 'r') as infile:
        sequences = []
        current_seq = []
        header = ""
        
        for line in infile:
            if line.startswith('>'):
                if current_seq:
                    sequences.append((header, "".join(current_seq)))
                header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_seq:
            sequences.append((header, "".join(current_seq)))
        
        total_sequences = len(sequences)
        split_size = (total_sequences // num_files) + (1 if total_sequences % num_files else 0)
        
        for file_index in range(num_files):
            output_fasta = f"{output_prefix}_{file_index+1}.fasta"
            with open(output_fasta, 'w') as outfile:
                for header, sequence in sequences[file_index * split_size:(file_index + 1) * split_size]:
                    outfile.write(f"{header}\n{sequence}\n")

# Example usage
split_fasta('/Users/zainamoussa/Desktop/Leishmania_amazonensis_pipeline/BLAST_amazonensis_guides.fa', 'split_output')
