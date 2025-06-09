#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:04:03 2025

@author: zainamoussa
"""

from Bio import Entrez
from Bio import SeqIO
import csv

# Set your email to comply with NCBI's policy
Entrez.email = "zaina_moussa@berkeley.edu"

# Search NCBI for Leishmania viannia kinetoplast minicircle sequences
def search_ncbi():
    query = 'leishmania major kinetoplast AND 400:900[Sequence Length]'
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=5000)  # Adjust retmax as needed
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Fetch metadata for the given sequence IDs
def fetch_metadata(id_list):
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "gb"))
    handle.close()

    # Filter out records containing 'amazonensis', 'partial', or 'peruviana' in the description
    filtered_records = []
    filtered_ids = []  # Store the IDs of the valid sequences
    for record in records:
        description = record.description.lower()  # Convert to lowercase for case-insensitive matching
        if any(exclude_term in description for exclude_term in [ "mexicana","partial", "peruviana", "probe", "infantum", "chagasi","tropica"]):
            print(f"Excluding sequence: {record.id} - {record.description}")
        else:
            filtered_records.append(record)
            filtered_ids.append(record.id)
    
    return filtered_records, filtered_ids


# Fetch FASTA sequences for the given sequence IDs
def fetch_fasta(id_list, output_fasta):
    with open(output_fasta, "w") as fasta_file:
        handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        
        # Split by FASTA headers
        entries = fasta_data.split("\n>")
        formatted_entries = []
        
        for entry in entries:
            lines = entry.split("\n")
            header = lines[0]
            sequence = "\n".join(lines[1:])
            
            # Remove the "/1-715" part if present
            cleaned_header = header.split(" ")[0].split("/")[0]
            description = " ".join(header.split(" ")[1:])
            formatted_header = f"{cleaned_header} {description}".strip()
            
            formatted_entries.append(f">{formatted_header}\n{sequence}")
        
        fasta_file.write("\n\n".join(formatted_entries) + "\n\n")
        
    print(f"FASTA sequences saved to {output_fasta}")

    


# Extract strain information from the features
def get_strain(record):
    for feature in record.features:
        if feature.type == "source":
            strain = feature.qualifiers.get("strain", ["Unknown"])[0]
            return strain
    return "Unknown"



# Save metadata to a TSV file
def save_to_tsv(records, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # Write header
        writer.writerow(["accession", "description", "length", "feature types", "subtype", "strain", "topSubtype", "subtypeFolder", "spp", "species"])
        
        for record in records:
            feature_types = "; ".join(set(f.type for f in record.features))
            strain = get_strain(record)
            top_subtype = "TRUE"  # Column filled with "TRUE" as requested
            spp = "leishmania"
            species = "leishmania"
            subtype = " ".join(record.description.split()[:2])  # First two words of the description
            subtype_folder = "subtype_" + subtype
            writer.writerow([record.id, record.description, len(record.seq), feature_types, subtype, strain, top_subtype, subtype_folder, spp, species])
    print(f"Metadata saved to {output_file}")

def main():
    try:
        ids = search_ncbi()
        print(f"Found {len(ids)} sequences.")
        if ids:
            # Fetch metadata and save to TSV
            records, filtered_ids = fetch_metadata(ids)
            save_to_tsv(records, "Desktop/Leishmania_donovani_metadata.tsv")
            
            # Fetch FASTA sequences and save to a file using the filtered IDs
            fetch_fasta(filtered_ids, "Desktop/Leishmania_donovani_sequences.fasta")
        else:
            print("No sequences found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()