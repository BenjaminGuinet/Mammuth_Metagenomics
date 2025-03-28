import argparse
import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm  # Import tqdm

def create_sub_msa(input_file, species_to_remove, output_dirm, ali_dir):
    # Read the MSA file
    alignment = AlignIO.read(input_file, "fasta")
    # Create a set of species to remove for faster lookup
    species_to_remove_set = set(species_to_remove)
    # Check if the species to remove are in the alignment
    species_in_alignment = [record.id for record in alignment]
    species_present = [species for species in species_to_remove_set if species in species_in_alignment]
    # Skip the file if no species are present
    if not species_present:
        print(f"Skipping file {input_file} as none of the species are present.")
        return
    # Iterate through species in the alignment
    for species in species_present:
        # Filter the alignment to include only the species in question and others
        filtered_species = [record for record in alignment if species in record.id or not any(spec in record.id for spec in species_to_remove_set)]   
        # Find the sequence with Ns
        ref_seq_record = next(rec for rec in filtered_species if rec.id == species)
        ref_seq = ref_seq_record.seq
        # Identify non-N positions
        valid_positions = [i for i, base in enumerate(ref_seq) if base != "N"]
        # Filter each sequence based on valid positions
        filtered_species = [
            SeqRecord(Seq("".join(seq[i] for i in valid_positions)), id=rec.id, name=rec.name, description=rec.description, dbxrefs=rec.dbxrefs)
            for rec in filtered_species
            for seq in [rec.seq]]
        # Convert the filtered records into a MultipleSeqAlignment object
        filtered_alignment = MultipleSeqAlignment(filtered_species)
        for record in filtered_alignment:
            # Modify the ID by splitting at the first period and taking the first part
            record.id = record.id.split('.')[0]
            record.description = record.description.split('.')[0]
        # Prepare output filename based on species which is the name of the file changing ALL_cat by the name of the species
        output_filename  = os.path.join(output_dirm, os.path.basename(input_file).replace("ALL_cat", "Sub_"+species))
        # Write the filtered alignment to the new file
        with open(output_filename, "w") as out_file:
            AlignIO.write(filtered_alignment, out_file, "fasta")
        # print the size of the filtered alignment
        #print(f"Filtered alignment size: {filtered_alignment.get_alignment_length()}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Create sub-alignments based on species to remove")
    # Adding arguments for input file, species, and output directory
    parser.add_argument("-s", "--species", nargs='+', help="Species to remove (space-separated)", required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save the output sub-alignments", required=True)
    parser.add_argument("-ali_dir", "--alignment_dir", help="Directory of the alignments", required=True)
    args = parser.parse_args()
    # Get species to remove from the command-line arguments
    species_to_remove = args.species
    # Get the directory path for the input files
    ali_dir = args.alignment_dir
    # Find all files matching the pattern
    input_files = [os.path.join(ali_dir, f) for f in os.listdir(ali_dir) if f.startswith("ALL_cat_edit0-2") ]
    # Use tqdm to show progress for processing files
    for input_file in tqdm(input_files, desc="Processing files", unit="file", total=len(input_files)):
        create_sub_msa(input_file, species_to_remove, args.output_dir, args.alignment_dir)

if __name__ == "__main__":
    main()
