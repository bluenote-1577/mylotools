#!/usr/bin/env python3
import argparse
import os
import subprocess
import tempfile
from Bio import SeqIO
from io import StringIO
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from plotly.offline import plot

def extract_contig(contig_id, fasta_file):
    """Extract a specific contig from the FASTA file."""
    print(f"Extracting contig {contig_id} from {fasta_file}...")
    cmd = f"grep -A1 '{contig_id}' {fasta_file}"
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        output = result.stdout
        # Process the output to get the proper FASTA format
        if not output:
            print(f"Error: Contig {contig_id} not found in {fasta_file}")
            return None
            
        lines = output.strip().split('\n')
        if len(lines) < 2:
            print(f"Error: Invalid output format when extracting contig {contig_id}")
            return None
            
        header = lines[0]
        if not header.startswith('>'):
            header = '>' + header
            
        sequence = lines[1]
        return f"{header}\n{sequence}"
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(f"Error message: {e}")
        return None


def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    sequence = sequence.upper()
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_count = g_count + c_count
    gminusc_count = g_count - c_count
    if len(sequence) > 0:
        return ((gc_count / len(sequence)) * 100, (gminusc_count / len(sequence)) * 100)
    else:
        return (0,0)

def calculate_gc_content_windows(sequence, window_size=1000, step_size=None):
    """Calculate GC content across a genome using a sliding window approach."""
    if step_size is None:
        step_size = window_size // 2
        
    sequence = sequence.upper()
    seq_length = len(sequence)
    positions = []
    gc_contents = []
    gc_skew_contents = []
    
    for i in range(0, seq_length - window_size + 1, step_size):
        window = sequence[i:i+window_size]
        gc, gc_skew = calculate_gc_content(window)
        positions.append(i + window_size // 2)  # Middle position of the window
        gc_contents.append(gc)
        gc_skew_contents.append(gc_skew)
        
    return positions, gc_contents, gc_skew_contents

def parse_gfa_for_coverage(contig_id, gfa_file):
    """Extract coverage information from GFA file for a specific contig."""
    print(f"Extracting coverage data for contig {contig_id} from {gfa_file}...")
    cmd = f"grep '{contig_id}' {gfa_file}"
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        output = result.stdout
        
        if not output:
            print(f"Error: Contig {contig_id} not found in {gfa_file}")
            return None, None, None, None, None, None, None, None, None
            
        positions = []
        dp1_values = []
        dp2_values = []
        dp3_values = []
        read_ids = []
        read_lens = []
        ol_lens = []
        snp_shares = []
        snp_diffs = []
        current_pos = 0
        
        for line in output.strip().split('\n'):
            if not line.startswith('a'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 8:  # Ensure we have enough columns
                continue
            
            # Extract the read ID from column 5
            read_id = parts[4] if len(parts) > 4 else "unknown_read"
                
            # Extract DP values and additional metrics
            dp_part = parts[-1]
            try:
                # Parse metrics using a more flexible approach to handle the extended format
                metrics = {}
                for metric in dp_part.split(','):
                    if ':' in metric:
                        key, value = metric.split(':', 1)
                        try:
                            # Try to convert to integer if possible
                            metrics[key] = int(value)
                        except ValueError:
                            # If not an integer, keep as string
                            metrics[key] = value
                
                # Extract required values with defaults if not present
                dp1 = metrics.get('DP1', 0)
                dp2 = metrics.get('DP2', 0)
                dp3 = metrics.get('DP3', 0)
                dp1_values.append(dp1)
                dp2_values.append(dp2)
                dp3_values.append(dp3)
                read_ids.append(read_id)
                length = int(parts[6])
                current_pos += length

                read_len = metrics.get('READ_LEN', 0)
                ol_len = metrics.get('OL_LEN_NEXT', 0)
                snp_share = metrics.get('SNP_SHARE_NEXT', 0)
                snp_diff = metrics.get('SNP_DIFF_NEXT', 0)
                
                # Get length for position calculation
                
                positions.append(current_pos)
                read_lens.append(read_len)
                ol_lens.append(ol_len)
                snp_shares.append(snp_share)
                snp_diffs.append(snp_diff)
            except (IndexError, ValueError) as e:
                print(f"Error parsing line: {line}")
                print(f"Error details: {e}")
                continue
                
        return positions, dp1_values, dp2_values, dp3_values, read_ids, read_lens, ol_lens, snp_shares, snp_diffs
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(f"Error message: {e}")
        return None, None, None, None, None, None, None, None, None
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {cmd}")
        print(f"Error message: {e}")
        return None, None, None, None

def create_interactive_plot(contig_id, gc_data, coverage_data, output_file):
    """Create an interactive HTML plot with GC content and coverage data."""
    positions_gc, gc_contents, gc_skew_contents = gc_data

    positions_cov, dp1_values, dp2_values, dp3_values, read_ids, read_lens, ol_lens, snp_shares, snp_diffs = coverage_data
    
    # Create subplots: GC content on top, coverage in middle, edge info on bottom
    fig = make_subplots(rows=4, cols=1, 
                        shared_xaxes=True, 
                        vertical_spacing=0.1,
                        subplot_titles=("GC Content", "GC Skew", "Coverage", "Read Overlap Information"),
                        row_heights=[0.3, 0.3, 0.5, 0.5])
    
    # Add GC content trace
    fig.add_trace(
        go.Scatter(
            x=positions_gc, 
            y=gc_contents,
            mode='lines',
            name='GC Content',
            line=dict(color='blue', width=1)
        ),
        row=1, col=1
    )
    
    # Add a horizontal line for mean GC content
    mean_gc = sum(gc_contents) / len(gc_contents)
    fig.add_trace(
        go.Scatter(
            x=[min(positions_gc), max(positions_gc)],
            y=[mean_gc, mean_gc],
            mode='lines',
            name=f'Mean GC: {mean_gc:.2f}%',
            line=dict(color='red', width=1, dash='dash')
        ),
        row=1, col=1
    )

    cumulative_gc_skew = []
    last = None
    for gc_skew in gc_skew_contents:
        if last == None:
            skew_val = gc_skew
        else:
            skew_val = skew_val + gc_skew
        cumulative_gc_skew.append(skew_val)
        last = skew_val

    # Add GC content trace
    fig.add_trace(
        go.Scatter(
            x=positions_gc, 
            y=cumulative_gc_skew,
            mode='lines',
            name='GC Skew',
            line=dict(color='blue', width=1)
        ),
        row=2, col=1
    )

    
    # Add coverage traces (only if we have coverage data)
    if positions_cov and dp1_values:
        # DP1 with lines and markers
        fig.add_trace(
            go.Scatter(
                x=positions_cov, 
                y=dp1_values,
                mode='lines+markers',
                name='DP1',
                line=dict(color='green'),
                marker=dict(size=8, symbol='circle', color='green'),
                hoverinfo='text',
                hovertext=[f'Position: {pos}<br>DP1: {dp}<br>Read ID: {read_id}' 
                           for pos, dp, read_id in zip(positions_cov, dp1_values, read_ids)]
            ),
            row=3, col=1
        )
        
        # DP2 with lines and markers
        fig.add_trace(
            go.Scatter(
                x=positions_cov, 
                y=dp2_values,
                mode='lines+markers',
                name='DP2',
                line=dict(color='orange'),
                marker=dict(size=8, symbol='circle', color='orange'),
                hoverinfo='text',
                hovertext=[f'Position: {pos}<br>DP2: {dp}<br>Read ID: {read_id}' 
                           for pos, dp, read_id in zip(positions_cov, dp2_values, read_ids)]
            ),
            row=3, col=1
        )
        
        # DP3 with lines and markers
        fig.add_trace(
            go.Scatter(
                x=positions_cov, 
                y=dp3_values,
                mode='lines+markers',
                name='DP3',
                line=dict(color='purple'),
                marker=dict(size=8, symbol='circle', color='purple'),
                hoverinfo='text',
                hovertext=[f'Position: {pos}<br>DP3: {dp}<br>Read ID: {read_id}' 
                           for pos, dp, read_id in zip(positions_cov, dp3_values, read_ids)]
            ),
            row=3, col=1
        )
        
        # Add edge overlap information scatter plot (excluding the last node)
        if len(positions_cov) > 1:  # Only if we have more than one position
            # Create a color scale based on SNP differences
            max_snp_diff = max(snp_diffs[:-1]) if snp_diffs[:-1] else 1  # Avoid division by zero
            colors = []
            for snp_diff in snp_diffs[:-1]:
                # Generate color from green (0) to red (max diff)
                # Lower SNP differences are better (more green)
                normalized_diff = snp_diff / max_snp_diff if max_snp_diff > 0 else 0
                colors.append(f'rgb({int(255 * normalized_diff)}, {int(255 * (1 - normalized_diff))}, 0)')
            
            fig.add_trace(
                go.Scatter(
                    x=positions_cov[:-1],  # Exclude the last node
                    y=ol_lens[:-1],        # Exclude the last node
                    mode='lines+markers',
                    name='Read Overlaps',
                    line=dict(color='black'),
                    marker=dict(
                        size=12,
                        color=colors,
                        symbol='circle',
                        line=dict(width=1, color='black')
                    ),
                    hoverinfo='text',
                    hovertext=[f'Position: {pos}<br>Read ID: {read_id}<br>Read Length: {read_len}<br>'
                               f'Overlap Length: {ol_len}<br>Shared SNPs: {snp_share}<br>'
                               f'Differing SNPs: {snp_diff}'
                               for pos, read_id, read_len, ol_len, snp_share, snp_diff 
                               in zip(positions_cov[:-1], read_ids[:-1], read_lens[:-1], 
                                     ol_lens[:-1], snp_shares[:-1], snp_diffs[:-1])]
                ),
                row=4, col=1
            )
            
            # Add a color bar for SNP differences
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(
                        colorscale=[
                            [0, 'green'],
                            [1, 'red']
                        ],
                        showscale=True,
                        colorbar=dict(
                            title='SNP Differences',
                            thickness=15,
                            len=0.5,
                            y=0.15,
                            yanchor='middle'
                        ),
                        cmin=0,
                        cmax=max_snp_diff
                    ),
                    showlegend=False,
                ),
                row=4, col=1
            )
    
    # Update layout
    fig.update_layout(
        title=f'Genome Analysis for Contig: {contig_id}',
        height=1000,  # Increased height to accommodate the third plot
        width=1200,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Update x-axis labels
    fig.update_xaxes(title_text="Position (bp)", row=4, col=1)
    
    # Update y-axis labels
    fig.update_yaxes(title_text="GC Content (%)", row=1, col=1)
    fig.update_yaxes(title_text="GC Skew", row=2, col=1)
    fig.update_yaxes(title_text="Coverage", row=3, col=1)
    fig.update_yaxes(title_text="Overlap Length (bp)", row=4, col=1)
    
    # Save as an HTML file
    plot(fig, filename=output_file, auto_open=False)
    print(f"Interactive plot saved as {output_file}")
    
    return fig
    
    # Update layout
    fig.update_layout(
        title=f'Genome Analysis for Contig: {contig_id}',
        height=800,
        width=1200,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Update x-axis labels
    fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
    
    # Update y-axis labels
    fig.update_yaxes(title_text="GC Content (%)", row=1, col=1)
    fig.update_yaxes(title_text="Coverage", row=2, col=1)
    
    # Save as an HTML file
    plot(fig, filename=output_file, auto_open=False)
    print(f"Interactive plot saved as {output_file}")
    
    return fig

def main(args):
        # Set default output file name if not specified
    for contig_id in args.contig_ids:
        if args.output is None:
            args.output = f"{contig_id}_analysis.html"
        
        # Extract the contig sequence
        contig_fasta = extract_contig(contig_id, args.fasta)
        if not contig_fasta:
            return
        
        # Parse the contig sequence
        record = list(SeqIO.parse(StringIO(contig_fasta), "fasta"))[0]
        sequence = str(record.seq)
        
        print(f"Analyzing contig {contig_id} (length: {len(sequence)} bp)")
        
        # Calculate GC content
        gc_data = calculate_gc_content_windows(sequence, args.window_size, args.step_size)
        
        # Get coverage data
        coverage_data = parse_gfa_for_coverage(contig_id, args.gfa)
        
        # Create the interactive plot
        create_interactive_plot(contig_id, gc_data, coverage_data, args.output)
    
        print(f"Analysis complete. Open {args.output} in a web browser to view the interactive plot.")
