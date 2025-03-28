import argparse
from Bio import SeqIO

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Remove sequences from a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("-delete", nargs='+', help="List of sequence IDs to remove", required=True)
    parser.add_argument("-o", "--output", help="Output FASTA file", required=True)

    # Parse arguments
    args = parser.parse_args()

    # Load the FASTA file
    sequences = list(SeqIO.parse(args.input_fasta, "fasta"))

    # Filter out sequences with IDs in the delete list
    filtered_sequences = [seq for seq in sequences if seq.id not in args.delete]

    # Write the filtered sequences to the output file
    SeqIO.write(filtered_sequences, args.output, "fasta")

    print(f"Sequences removed and saved to {args.output}")

if __name__ == "__main__":
    main()

