import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
import io
import time

# Define exclusion keywords
EXCLUDE_KEYWORDS = {
    "Leishmania viannia",
    "Leishmania braziliensis",
    "Leishmania panamensis",
    "Leishmania guyanensis",
    "Leishmania lainsoni"
}

# Load DNA sequences from a text file
def load_sequences(txt_file, column_name="target"):
    # Read the text file into a DataFrame
    df = pd.read_csv(txt_file, sep='\t')  # Assuming tab-separated values
    return df[column_name].dropna().tolist()

# Run BLAST for each sequence
def run_blast(sequence):
    try:
        print(f"Running BLAST for sequence: {sequence[:30]}...")  # Show part of the sequence
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        return result_handle.read()
    except Exception as e:
        print(f"BLAST failed for sequence: {sequence[:30]}... Error: {str(e)}")
        return None

# Parse and filter BLAST results
def parse_blast_results(blast_xml):
    if not blast_xml:
        return []
    
    blast_records = NCBIXML.read(io.StringIO(blast_xml))
    filtered_results = []

    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            identity = (hsp.identities / hsp.align_length) * 100  # Calculate % identity
            query_coverage = (hsp.query_end - hsp.query_start + 1) / len(hsp.query) * 100  # Query coverage
            
            if identity >=90 and query_coverage >=90 and not any(keyword in alignment.title for keyword in EXCLUDE_KEYWORDS):
                filtered_results.append({
                    "Title": alignment.title,
                    "Length": alignment.length,
                    "Score": hsp.score,
                    "E-value": hsp.expect,
                    "Identity (%)": identity,
                    "Query Coverage (%)": query_coverage
                })
    
    return filtered_results

# Save results incrementally
def save_results(all_results, output_file="blast_results.xlsx"):
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_excel(output_file, index=False)
        print(f"Results saved to {output_file}")
    else:
        print("No valid results to save.")

# Main function with batching
def main(txt_file, output_file="blast_results.xlsx", batch_size=10):
    sequences = load_sequences(txt_file)
    print(sequences[0])
    all_results = []

    # Run in batches of `batch_size`
    '''
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}...")
        
        for seq in batch:
            blast_xml = run_blast(seq)
            filtered_hits = parse_blast_results(blast_xml)
            all_results.extend(filtered_hits)
        
        # Save progress after each batch
        save_results(all_results, output_file)

        # Pause to prevent too many requests in a short period
        print("Waiting for a few seconds before the next batch...")
        time.sleep(5)  # Adjust time (in seconds) if necessary
        '''
    blast_xml = run_blast(sequences[0])
    all_results.extend
# Example usage
if __name__ == "__main__":
    # Specify the full path to your text file here
    input_file_path = "/Users/zainamoussa/Desktop/Leishmania_pipeline/Cas12_cross_reactivity_guide_results_summary.txt"  # Update with your actual path
    main(input_file_path)
