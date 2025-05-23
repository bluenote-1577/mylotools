import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import gzip
import re

def parse_overlaps_file(overlaps_file):
    """Parse the overlaps file to extract read connections."""
    overlaps = []
    
    # Check if the file is gzipped
    open_func = gzip.open if overlaps_file.endswith('.gz') else open
    mode = 'rt' if overlaps_file.endswith('.gz') else 'r'
    
    with open_func(overlaps_file, mode) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 10:  # Basic validation
                continue
                
            read1_id = parts[2]  # DRR582205.796678
            read2_id = parts[3]  # DRR582205.1433003
            
            # Extract FSV (overlap quality)
            fsv_match = re.search(r'fsv:(\d+\.?\d*)', line)
            fsv = float(fsv_match.group(1)) if fsv_match else 0
            
            # Extract SHARE (overlap size)
            ol_match = re.search(r'(\d+-\d+)', parts[10])
            ol = ol_match.group(1).split('-')
            ol_len = int(ol[1]) - int(ol[0])
            
            overlaps.append({
                'read1': read1_id,
                'read2': read2_id,
                'fsv': fsv,
                'ol_len': ol_len
            })
    
    return overlaps

def parse_gfa_file(gfa_file):
    """Parse GFA file to extract contigs and their constituent reads."""
    contigs = {}
    current_contig = None
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):  # Segment line
                parts = line.strip().split('\t')
                contig_id = parts[1]
                contigs[contig_id] = []
                current_contig = contig_id
            elif line.startswith('a') and current_contig:  # Alignment line in some GFA variants
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    read_info = parts[4]  # This field should contain the read ID
                    read_id = read_info.split()[0]  # Extract just the ID part
                    contigs[current_contig].append(read_id)
    
    return contigs

def filter_contigs(contigs, contig_ids=None, max_reads=None):
    """Filter contigs by ID and/or limit the number of reads per contig."""
    filtered_contigs = {}
    
    # Filter by contig IDs if specified
    if contig_ids:
        for contig_id in contig_ids:
            if contig_id in contigs:
                filtered_contigs[contig_id] = contigs[contig_id]
    else:
        filtered_contigs = contigs
    
    # Limit number of reads per contig if specified
    if max_reads:
        for contig_id in filtered_contigs:
            filtered_contigs[contig_id] = filtered_contigs[contig_id][:max_reads]
    
    return filtered_contigs

def plot_contig_overlaps(contigs, overlaps, max_reads_per_contig, output_file=None):
    """Create a visualization of contig overlaps."""
    # Limit the number of reads per contig for visualization
    filtered_contigs = {}
    for contig_id, reads in contigs.items():
        filtered_contigs[contig_id] = reads[:max_reads_per_contig]
    
    # Create a map of read IDs to their positions (contig and index)
    read_positions = {}
    for contig_idx, (contig_id, reads) in enumerate(filtered_contigs.items()):
        for read_idx, read_id in enumerate(reads):
            read_positions[read_id] = (contig_idx, read_idx)
    
    # Set up the plot
    reads_per_contig = [len(reads) for _, reads in contigs.items()]
    max_reads = max(reads_per_contig)
    fig, ax = plt.subplots(figsize=(min(max(12, max_reads/10), 24), len(filtered_contigs) * 2))
    
    # Draw reads as boxes
    box_width = 0.8
    box_height = 0.4
    
    # Plot each contig as a row
    for contig_idx, (contig_id, reads) in enumerate(filtered_contigs.items()):
        y_pos = contig_idx * 2  # Vertical position for this contig
        
        # Draw the contig label
        ax.text(-1, y_pos, contig_id, va='center', ha='right', fontsize=10)
        
        # Draw each read as a box
        for read_idx, read_id in enumerate(reads):
            x_pos = read_idx
            
            # Create a box for the read
            rect = patches.Rectangle((x_pos, y_pos - box_height/2), box_width, box_height, 
                                     linewidth=1, edgecolor='black', facecolor='lightgray')
            ax.add_patch(rect)
            
            # Add read ID label (optional, might be too cluttered)
            if len(reads) < 10:  # Only show labels if not too many reads
                ax.text(x_pos + box_width/2, y_pos, read_id.split('.')[-1], 
                        va='center', ha='center', fontsize=8, rotation=90)
    
    # Draw connections between reads based on overlaps
    for overlap in overlaps:
        read1 = overlap['read1']
        read2 = overlap['read2']
        
        # Skip if either read is not in our visualization
        if read1 not in read_positions or read2 not in read_positions:
            continue
        
        contig1_idx, read1_idx = read_positions[read1]
        contig2_idx, read2_idx = read_positions[read2]
        
        # Calculate positions
        x1 = read1_idx + box_width/2
        y1 = contig1_idx * 2
        x2 = read2_idx + box_width/2
        y2 = contig2_idx * 2
        
        ol_len = overlap['ol_len']
        
        # Normalize share for line thickness (capped at 5 for visibility)
        thickness = min(2, ol_len / 20000)

        fsv = overlap['fsv']

        # Color based on FSV (red for low, green for high)
        if fsv < 98:
            color = 'purple'
            alpha = 0.3
        elif fsv < 99.8:
            color = 'orange'
            alpha = 0.5
        elif fsv < 99.99:
            color = 'blue'
            alpha = 0.7
        else:
            color = 'green'
            alpha = 0.9
        
        # Draw the connection
        if contig1_idx == contig2_idx:
            # Connection within the same contig (curved line above)
            ax.plot([x1, (x1+x2)/2, x2], 
                    [y1, y1 + 0.5, y1], 
                    color=color, alpha=alpha, linewidth=thickness)
        else:
            # Connection between different contigs (straight line)
            ax.plot([x1, x2], [y1, y2], color=color, alpha=alpha, 
                    linewidth=thickness, linestyle='--')
    
    # Set plot limits and remove axes
    ax.set_xlim(-1.5, max([len(reads) for reads in filtered_contigs.values()]))
    ax.set_ylim(-1, len(filtered_contigs) * 2)
    
    # Add legend for FSV colors
    legend_elements = [
        patches.Patch(facecolor='green', alpha=0.9, label='FSV = 100'),
        patches.Patch(facecolor='blue', alpha=0.7, label='99.8 ≤ FSV < 99.99'),
        patches.Patch(facecolor='orange', alpha=0.5, label='98 ≤ FSV < 99.8'),
        patches.Patch(facecolor='purple', alpha=0.3, label='FSV < 98')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add title and labels
    plt.title('Contig Read Overlaps Visualization')
    plt.xlabel('Read Position in Contig')
    plt.yticks([])
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Visualization saved to {output_file}")
    else:
        plt.show()

def main(args):

    # Parse input files
    overlaps = parse_overlaps_file(args.overlaps)
    contigs = parse_gfa_file(args.gfa)
    
    # Filter contigs if specified
    contig_ids = args.contigs
    filtered_contigs = filter_contigs(contigs, contig_ids, args.max_reads)
    
    # Create visualization
    if args.output == "" or args.output == None:
        print("Outputting visualization to terminal. Specify --output to plot to file.")

    plot_contig_overlaps(filtered_contigs, overlaps, args.max_reads, args.output)
