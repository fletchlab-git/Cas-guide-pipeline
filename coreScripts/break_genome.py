import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default="ref_data/NC_045512v2.fa",
                    help="genome .fa file name")
parser.add_argument("-d", "--scriptdir", default=".",
                    help="dummy argument for compatibility with break_genomes.py")
parser.add_argument("-w", "--window", type=int, default=20,
                    help="window size")
parser.add_argument("-e", "--enzyme", type=str, help="Cas enzyme type")  # "Cas13a" or "Cas12"
parser.add_argument("-p", "--pfs_length", type=int, default=4,
                    help="length of PFS/PAM to evaluate")
parser.add_argument("-g", "--genome_strand", type=str, default="+",
                    help="strandedness of viral genome")
parser.add_argument("-s", "--strand", type=str, default="+",
                    help="strand to target")
parser.add_argument("-o", "--out", type=str, default="break_genome_output",
                    help="output directory")
args = parser.parse_args()

if args.enzyme is None:
    print("\nNO ENZYME SPECIFIED\n")
    exit(1)
if not os.path.exists(args.out):
    os.makedirs(args.out)


genome_seq = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

windows = []
for seq_id, seq in genome_seq.items():
    seq_str = str(seq.seq)
    num_windows = len(seq_str) - args.window - 3
    for i in range(num_windows):
        windows.append({
            'segment': seq_id.split(" ")[0],
            'start': i + 1,
            'target': seq_str[i: i + args.window]
        })

df = pd.DataFrame(windows)
num_rows, num_columns = df.shape
if num_rows == 0:
    print(f"No windows found in genome {args.input}. Exiting...")
    exit(1)
if args.enzyme == "Cas13a":
    if args.strand == args.genome_strand:
        df['spacer'] = df['target'].apply(lambda x: Seq(x).reverse_complement().transcribe())
        df['strand'] = args.genome_strand
    else:
        df['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
        df['target'] = df['target'].apply(lambda x: Seq(x).reverse_complement())
        df['strand'] = args.strand

elif args.enzyme == "Cas12":
    # Cas12 targets DNA
    df['spacer'] = df['target'].apply(lambda x: Seq(x).reverse_complement().transcribe())
    df['strand'] = args.genome_strand

    df_minusStrand = df.copy()
    df_minusStrand['target'] = df['target'].apply(lambda x: Seq(x).reverse_complement()) #fix this!!!
    df_minusStrand['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
    df_minusStrand['strand'] = '-'
    df = pd.concat([df, df_minusStrand], ignore_index=True)

# Apply filters based on tag/antitag and/or PAM if needed (Cas13a or Cas12)
if args.enzyme == "Cas13a":
    # Cas13a-specific tag/antitag processing
    df['antitag'] = df['start'].apply(lambda x: genome_seq[seq_id].seq[x + args.window: x + args.window + 3].replace("T", "U"))

if args.enzyme == "Cas12":
    # Cas12 targets DNA, initialize PAM column with None
   # Process PAM sequences for Cas12 (4-nucleotide PAM)
    df['PAM'] = None

    # Process PAM for the positive strand ("+")
    for i, row in df[df['strand'] == "+"].iterrows():
        x = row['start']
        segment_name = row['segment']
        genome_sequence = genome_seq[segment_name].seq  # Get the genome sequence for the segment
        if x + args.window + 4 <= len(genome_sequence):  # Ensure we're within bounds for 4 nucleotides
            PAM = genome_sequence[x + args.window-1: x + args.window + 4-1]  # Extract 4-nt PAM sequence
            df.at[i, 'PAM'] = str(Seq(PAM).reverse_complement())  # Reverse complement the PAM
        elif x+args.window+4 >= len(genome_sequence): #not gonna worry about circular DNA for now smhh too much
            loop = (x+args.window+4)-len(genome_sequence) #this is circular dna so go around the other side
            PAM_end = genome_sequence[x + args.window-1:]
            PAM_loop = genome_sequence[0:loop]
        else:
            df.at[i, 'PAM'] = None  # Set to None if out of bounds

    # Process PAM for the negative strand ("-")
    for i, row in df[df['strand'] == "-"].iterrows():
        x = row['start']
        segment_name = row['segment']
        genome_sequence = genome_seq[segment_name].seq  # Get the genome sequence for the segment
        if x - 4 >= 0:  # Ensure we're within bounds for the negative strand
            PAM = genome_sequence[x -4-1: x - 1]  # Extract PAM for the negative strand (4 nucleotides)
            df.at[i, 'PAM'] = str(PAM)  # Set the PAM for the negative strand
          
    # Filter PAM sequences: length >= 4 and start with "TTT"
    df = df[df['PAM'].str.len() >= 4]
    df = df[df['PAM'].str.startswith("TTT")]


# Assuming 'windows' is a pandas DataFrame with the equivalent structure in Python.
# Also assuming you have the 'genome_seq' as a dictionary where keys are segment names 
# and values are Seq objects from BioPython.



df['GC_content'] = df['spacer'].apply(lambda x: round(gc_fraction(Seq(x)),2))
df['A_content'] = df['spacer'].apply(lambda x: x.upper().count('A')/args.window)

df.to_csv(f"{args.out}/windows.txt", sep='\t', index=False)

df['target'].to_csv(f"{args.out}/targets.txt", index=False)
df['spacer'].to_csv(f"{args.out}/spacers.txt", index=False)

with open(f"{args.out}/targets.fa", 'w') as f:
    for idx, row in df.iterrows():
        f.write(f">{row['segment']}_{row['start']}\n")
        f.write(f"{row['target']}\n")
