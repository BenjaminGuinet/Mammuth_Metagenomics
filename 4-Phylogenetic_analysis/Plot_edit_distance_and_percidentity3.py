import os
import subprocess
import pysam
import matplotlib.pyplot as plt
import statistics
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages

def get_identity_distribution_and_avg_mapq(bam_file):
    identity_values = []
    edit_distances = []
    filtered_identity_values = []
    filtered_edit_distances = []
    total_nb_reads = 0
    total_mapq = []
    filtered_mapq = []
    
    if not os.path.exists(bam_file + ".bai"):
        command = ["samtools index", bam_file]
        subprocess.run(" ".join(command), shell=True, check=True)
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                edit_distance = read.get_tag("NM")  # NM tag represents edit distance
                read_length = read.query_length
                percentage_identity = ((read_length - edit_distance) / read_length) * 100
                identity_values.append(percentage_identity)
                edit_distances.append(edit_distance)
                total_mapq.append(read.mapping_quality)
                total_nb_reads += 1
                
                # Filtered lists for edit distance <= 2
                if edit_distance <= 2:
                    filtered_identity_values.append(percentage_identity)
                    filtered_edit_distances.append(edit_distance)
                    filtered_mapq.append(read.mapping_quality)
    
    avg_mapq = sum(total_mapq) / len(total_mapq) if total_mapq else 0
    avg_identity = sum(identity_values) / len(identity_values) if identity_values else 0
    filtered_avg_identity = sum(filtered_identity_values) / len(filtered_identity_values) if filtered_identity_values else 0
    filtered_avg_mapq = sum(filtered_mapq) / len(filtered_mapq) if filtered_mapq else 0
    
    return identity_values, edit_distances, filtered_identity_values, filtered_edit_distances, total_nb_reads, avg_mapq, avg_identity, filtered_avg_identity, filtered_avg_mapq

def process_bam_file(bam_file, pdf_pages):
    identity_values, edit_distances, filtered_identity_values, filtered_edit_distances, total_nb_reads, avg_mapq, avg_identity, filtered_avg_identity, filtered_avg_mapq = get_identity_distribution_and_avg_mapq(bam_file)
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))  # Single row, four columns
    plt.subplots_adjust(wspace=0.5)
    
    max_edit_distance = max(edit_distances) if edit_distances else 10
    
    # Edit Distance Histogram (Before Filtering)
    axes[0].hist(edit_distances, bins=50, color='blue', alpha=0.7)
    axes[0].set_title("Edit Distance (Before)")
    axes[0].set_xlabel("Edit Distance (NM tag)")
    axes[0].set_ylabel("Count")
    axes[0].set_xlim(0, 10)
    axes[0].text(0.95, 0.95, f'Avg MAPQ: {avg_mapq:.2f}', transform=axes[0].transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
    
    # Percentage Identity Histogram (Before Filtering)
    axes[1].hist(identity_values, bins=50, color='green', alpha=0.7)
    axes[1].set_title("% Identity (Before)")
    axes[1].set_xlim(85, 100)
    axes[1].set_xlabel("% Identity")
    axes[1].set_ylabel("Count")
    axes[1].text(0.95, 0.95, f'Avg Identity: {avg_identity:.2f}%', transform=axes[1].transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
    
    # Edit Distance Histogram (After Filtering: NM ≤ 2)
    axes[2].hist(filtered_edit_distances, bins=10, color='blue', alpha=0.7)
    axes[2].set_title("Edit Distance (After NM ≤ 2)")
    axes[2].set_xlabel("Edit Distance (NM tag)")
    axes[2].set_ylabel("Count")
    axes[2].set_xlim(0, 10)
    axes[2].text(0.95, 0.95, f'Avg MAPQ: {filtered_avg_mapq:.2f}', transform=axes[2].transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
    
    # Percentage Identity Histogram (After Filtering: NM ≤ 2)
    axes[3].hist(filtered_identity_values, bins=50, color='green', alpha=0.7)
    axes[3].set_title("% Identity (After NM ≤ 2)")
    axes[3].set_xlabel("% Identity")
    axes[3].set_ylabel("Count")
    axes[3].set_xlim(85, 100)
    axes[3].text(0.95, 0.95, f'Avg Identity: {filtered_avg_identity:.2f}%', transform=axes[3].transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))
    
    plt.suptitle(f"{bam_file}\n             Nb Reads Before: {total_nb_reads}                                                                                        Nb Reads filtred for phylogeny: {len(filtered_identity_values)}", weight='bold', fontsize=12)
    plt.tight_layout()
    pdf_pages.savefig(fig)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Process multiple BAM files and plot identity distributions.")
    parser.add_argument("--bam_files", nargs='+', type=str, required=True, help="List of BAM files.")
    parser.add_argument("--out", required=True, help="The pdf output file name")
    args = parser.parse_args()
    
    with PdfPages(args.out) as pdf_pages:
        for bam_file in tqdm(args.bam_files, desc="Processing BAM files"):
            process_bam_file(bam_file, pdf_pages)
    
if __name__ == "__main__":
    main()
