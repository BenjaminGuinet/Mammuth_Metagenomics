import argparse
from Bio import AlignIO
from tqdm import tqdm

def calculate_pairwise_identity(seq1, seq2):
    matches, valid_sites = 0, 0
    for a, b in zip(seq1, seq2):
        if a in ('-', 'N') or b in ('-', 'N'):
            continue
        valid_sites += 1
        if a == b:
            matches += 1
    return (matches / valid_sites * 100) if valid_sites > 0 else 0

def compute_identity_matrix(msa_file, output_file):
    alignment = AlignIO.read(msa_file, "fasta")
    sequences = [record.seq for record in alignment]
    seq_ids = [record.id for record in alignment]
    
    total_comparisons = (len(sequences) * (len(sequences) - 1)) // 2
    
    with open(output_file, "w") as out:
        out.write("Sequence1;Sequence2;Identity%\n")
        for i in tqdm(range(len(sequences)), desc="Computing identities"):
            for j in range(i + 1, len(sequences)):
                identity = calculate_pairwise_identity(sequences[i], sequences[j])
                out.write(f"{seq_ids[i]};{seq_ids[j]};{identity:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute pairwise sequence identity ignoring gaps and Ns.")
    parser.add_argument("-msa", required=True, help="Input multiple sequence alignment in FASTA format")
    parser.add_argument("-out", required=True, help="Output file for identity table")
    args = parser.parse_args()
    compute_identity_matrix(args.msa, args.out)
