from Bio import SeqIO
from pathlib import Path
import re

def extract_long_contigs(fasta_file, min_length, output_folder):
    """
    Extract contigs longer than min_length from a FASTA file and write each to separate files.
    Uses BioPython for robust FASTA parsing.
    
    Args:
        fasta_file (str): Path to input FASTA file
        min_length (int): Minimum contig length threshold
        output_folder (str): Directory to write individual contig files
    
    Returns:
        int: Number of contigs extracted
    """
    # Create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    contigs_extracted = 0
    
    # Parse FASTA file using BioPython
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > min_length:
            write_contig(record, output_folder)
            contigs_extracted += 1
    
    return contigs_extracted

def write_contig(record, output_folder):
    """
    Write a single contig to a file in the output folder.
    
    Args:
        record (Bio.SeqRecord): BioPython SeqRecord object
        output_folder (str): Output directory
    """
    # Clean record ID for filename (remove problematic characters)
    #safe_filename = re.sub(r'[^\w\s-]', '', record.id)
    #safe_filename = re.sub(r'[-\s]+', '_', safe_filename).strip('_')
    safe_filename = record.id
    
    # Create output file path
    output_path = Path(output_folder) / f"{safe_filename}.fa"
    
    # Write using BioPython's SeqIO
    SeqIO.write(record, output_path, "fasta")

# Alternative function that preserves original headers
def extract_long_contigs_preserve_headers(fasta_file, min_length, output_folder):
    """
    Same as extract_long_contigs but preserves full headers (description + id).
    
    Args:
        fasta_file (str): Path to input FASTA file
        min_length (int): Minimum contig length threshold
        output_folder (str): Directory to write individual contig files
    
    Returns:
        int: Number of contigs extracted
    """
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    contigs_extracted = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > min_length:
            # Use full description if available, otherwise use ID
            header = record.description if record.description else record.id
            safe_filename = re.sub(r'[^\w\s-]', '', header[:50])  # Limit length
            safe_filename = re.sub(r'[-\s]+', '_', safe_filename).strip('_')
            
            output_path = Path(output_folder) / f"{safe_filename}.fa"
            SeqIO.write(record, output_path, "fasta")
            contigs_extracted += 1
    
    return contigs_extracted

def main(args):
    num_extracted = extract_long_contigs(args.fasta, args.min_contig_length, args.output_folder)
    print(f"Extracted {num_extracted} contigs longer than {args.min_contig_length} bp to {args.output_folder}")
