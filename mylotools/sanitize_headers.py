#!/usr/bin/env python3
"""
Sanitize FASTA headers by replacing underscores with spaces (except the first one).
Creates a backup of the original file before modifying.
"""
from Bio import SeqIO
from pathlib import Path
import shutil


def sanitize_fasta_headers(fasta_file, backup=True):
    """
    Sanitize FASTA headers by replacing underscores with spaces.
    
    The contig ID (before first underscore) is preserved, but subsequent
    underscores are replaced with spaces for better readability.
    
    Example:
        u160616ctg_len-4164_circular-yes_depth-252-252-252_duplicated-no
        becomes:
        u160616ctg len-4164 circular-yes depth-252-252-252 duplicated-no
    
    Args:
        fasta_file (str): Path to input FASTA file
        backup (bool): Whether to create a backup (default: True)
    
    Returns:
        int: Number of headers modified
    """
    fasta_path = Path(fasta_file)
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"File not found: {fasta_file}")
    
    # Create backup if requested
    if backup:
        backup_path = fasta_path.with_suffix(fasta_path.suffix + '.bak')
        shutil.copy2(fasta_file, backup_path)
        print(f"Created backup: {backup_path}")
    
    # Read all records and modify headers
    modified_count = 0
    records = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        original_id = record.id
        original_desc = record.description
        
        # Split on first underscore only
        if '_' in original_id:
            parts = original_id.split('_', 1)
            base_id = parts[0]
            rest = parts[1] if len(parts) > 1 else ""
            
            # Replace remaining underscores with spaces
            new_rest = rest.replace('_', ' ')
            new_id = f"{base_id} {new_rest}" if new_rest else base_id
            
            if new_id != original_id:
                # Preserve any extra info after the ID (like "mult=...")
                # The description includes both ID and extra info
                extra_info = original_desc[len(original_id):].strip()
                record.id = new_id
                # Set description to just the extra info - SeqIO will prepend the ID
                record.description = extra_info
                modified_count += 1
                print(f"  {original_id} -> {new_id}")
        
        records.append(record)
    
    # Write modified records back to the same file
    SeqIO.write(records, fasta_file, "fasta")
    
    return modified_count


def main(args):
    """
    Main function for sanitize-headers command.
    
    Args:
        args: Parsed command-line arguments
    """
    print(f"\n{'='*60}")
    print(f"MyloTools FASTA Header Sanitizer")
    print(f"{'='*60}\n")
    
    print(f"Processing: {args.fasta}")
    
    try:
        modified = sanitize_fasta_headers(args.fasta, backup=not args.no_backup)
        
        print(f"\n{'='*60}")
        print(f"Sanitization complete!")
        print(f"{'='*60}")
        print(f"Modified {modified} header(s)")
        
        if not args.no_backup:
            print(f"Original backed up as: {args.fasta}.bak")
        print(f"Sanitized file: {args.fasta}\n")
        
    except Exception as e:
        print(f"\nError: {e}\n")
        return 1
    
    return 0
