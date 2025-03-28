import argparse
import pandas as pd
import re 
def parse_arguments():
    parser = argparse.ArgumentParser(description="Rename sequence IDs in ALN and Phylogeny files based on metadata.")
    parser.add_argument("metadata", help="Path to the metadata CSV file.")
    parser.add_argument("aln", help="Path to the ALN file.")
    parser.add_argument("phylogeny", help="Path to the phylogeny tree file.")
    return parser.parse_args()

def read_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path)
    metadata['Path'] = metadata['Path'].str.replace('./', '', regex=False).str.replace('_genomic.fna', '', regex=False)
    return metadata

def rename_sequences(metadata, aln_path, phylogeny_path):
    # Create a dictionary for quick lookup of new names
    name_mapping = {}
    name_counts = {}

    for _, row in metadata.iterrows():
        old_name = row['Path']
        new_name = row['Organism']

        if new_name in name_counts:
            name_counts[new_name] += 1
            new_name = f"{new_name}_X{name_counts[new_name]}"
        else:
            name_counts[new_name] = 0

        name_mapping[old_name] = new_name

    # Rename in ALN file
    with open(aln_path, 'r') as aln_file:
        aln_content = aln_file.readlines()

    with open(aln_path.replace('.aln', '_NewTipName.aln'), 'w') as new_aln_file:
        for line in aln_content:
            if line.startswith('>'):
                old_name = line[1:].strip()
                new_name = name_mapping.get(old_name, old_name)
                new_name=re.sub(" ","_",new_name)
                new_aln_file.write(f">{new_name}\n")
            else:
                new_aln_file.write(line)

    # Rename in Phylogeny file
    with open(phylogeny_path, 'r') as phylogeny_file:
        phylogeny_content = phylogeny_file.read()

    for old_name, new_name in name_mapping.items():
        phylogeny_content = phylogeny_content.replace(old_name, re.sub(" ","_",new_name))

    with open(phylogeny_path.replace('.treefile', '_NewTipName.treefile'), 'w') as new_phylogeny_file:
        new_phylogeny_file.write(phylogeny_content)

def main():
    args = parse_arguments()
    metadata = read_metadata(args.metadata)
    rename_sequences(metadata, args.aln, args.phylogeny)

if __name__ == "__main__":
    main()
