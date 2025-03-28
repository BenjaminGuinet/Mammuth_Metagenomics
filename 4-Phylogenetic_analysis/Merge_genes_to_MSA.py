import argparse
import os
from Bio import AlignIO

def merge_alignments(directory, output_file, partition_file,tag,  species_pair=None):
    alignments = []
    files = [file for file in os.listdir(directory) if file.startswith(tag) and file.endswith(".aln")]
    for file in files:
        alignment = AlignIO.read(os.path.join(directory, file), "fasta")
        for record in alignment:
            record.id = record.id.split(".")[0]
            record.description = record.description.split(".")[0]
        alignments.append(alignment)
    seq_ids = set()
    for alignment in alignments:
        seq_ids.update([record.id for record in alignment])
    merged_dict = {seq_id: [] for seq_id in seq_ids}
    partition_lines = []
    current_start = 1
    for file, alignment in zip(files, alignments):
        alignment_length = len(alignment[0])
        if species_pair:
            species_seqs = {sp: next((str(record.seq) for record in alignment if record.id == sp), None) for sp in species_pair}
            if None in species_seqs.values():
                print(f"Skipping {file} due to missing species.")
                continue
            valid_positions = [i for i in range(alignment_length) 
                               if all(species_seqs[sp][i] not in ('-', 'N') for sp in species_pair)]
        else:
            valid_positions = list(range(alignment_length))
        # Calculate length after filtering valid positions
        filtered_length = len(valid_positions)
        if filtered_length < 10:
            print(f"Skipping {file} as filtered length is less than 50bp.")
            continue  # Skip this file if the filtered length is less than 50bp
        for seq_id in seq_ids:
            seq = next((str(record.seq) for record in alignment if record.id == seq_id), None)
            if seq is None:
                seq = 'N' * alignment_length
            filtered_seq = ''.join(seq[i] for i in valid_positions)
            merged_dict[seq_id].append(filtered_seq)
        current_end = current_start + filtered_length - 1
        partition_lines.append(f"DNA, {file} = {current_start}-{current_end}")
        current_start = current_end + 1
    merged_alignment = [f">{seq_id}\n{''.join(sequences)}" for seq_id, sequences in merged_dict.items()]
    with open(output_file, "w") as output:
        output.write("\n".join(merged_alignment))
    with open(partition_file, "w") as partition:
        partition.write("\n".join(partition_lines))
    print(f"Merged alignment written to {output_file}")
    print(f"Partition file written to {partition_file}")

def main():
    parser = argparse.ArgumentParser(description="Merge sequence alignments into a single file and create partition file.")
    parser.add_argument("directory", help="Directory containing alignment files.")
    parser.add_argument("output_file", help="Path to the output merged alignment file.")
    parser.add_argument("partition_file", help="Path to the partition file.")
    parser.add_argument("--species_pair", help="List of two species for filtering.", nargs='+', default=None)
    parser.add_argument("--tag", help="Tag to find the begening of the ali file.")
    
    args = parser.parse_args()
    species_pair = args.species_pair
    merge_alignments(args.directory, args.output_file, args.partition_file, args.tag, species_pair)

if __name__ == "__main__":
    main()

"""
Eg usage:
python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Merge_genes_to_MSA.py ${output_main_dir}/5-align/Align-GENO3_1/ ${output_main_dir}/5-align/All_partitions_merged_filtred.aln ${output_main_dir}/5-align/All_partitions_merged_filtred.part
""" 
