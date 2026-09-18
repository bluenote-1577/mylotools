#!/usr/bin/env python3
"""
Report generation module for mylotools.
Creates comprehensive reports with individual contig plots and a summary visualization.
"""
import os
from pathlib import Path
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.offline import plot
import statistics
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from mylotools.plot import (
    calculate_gc_content_windows,
    create_interactive_plot
)


def extract_long_contigs(fasta_file, min_length, output_folder):
    """
    Extract contigs longer than min_length OR circular contigs from a FASTA file.
    
    Args:
        fasta_file (str): Path to input FASTA file
        min_length (int): Minimum contig length threshold
        output_folder (str): Directory to write individual contig files
    
    Returns:
        dict: Dictionary mapping contig_id -> {'length': int, 'filepath': str, 'full_id': str, 'is_circular': bool, 'metadata': dict}
    """
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    contig_info = {}
    
    # Parse FASTA file using BioPython
    for record in SeqIO.parse(fasta_file, "fasta"):
        full_contig_id = record.id
        contig_length = len(record.seq)
        
        # Parse metadata from contig ID
        # e.g., "u3601407ctg_len-2653912_circular-yes_depth-16-16-16_mult-1.00"
        metadata = {}
        is_circular = False
        
        if '_' in full_contig_id:
            parts = full_contig_id.split('_')
            for part in parts[1:]:  # Skip the base ID
                if '-' in part:
                    key, value = part.split('-', 1)
                    metadata[key] = value
                    if key == 'circular' and value == 'yes':
                        is_circular = True
        
        # Include contig if it's long enough OR if it's circular
        if len(record.seq) > min_length or is_circular:
            # Extract the base contig ID (before any underscore extensions)
            base_contig_id = full_contig_id.split('_')[0] if '_' in full_contig_id else full_contig_id
            
            # Write contig to individual file
            output_path = Path(output_folder) / f"{full_contig_id}.fa"
            SeqIO.write(record, output_path, "fasta")
            
            contig_info[base_contig_id] = {
                'length': contig_length,
                'filepath': str(output_path),
                'full_id': full_contig_id,
                'is_circular': is_circular,
                'metadata': metadata
            }
            
            circular_str = " [CIRCULAR]" if is_circular else ""
            print(f"Extracted contig {full_contig_id} ({contig_length:,} bp){circular_str}")
    
    return contig_info


def parse_gfa_for_all_contigs(gfa_file, contig_ids):
    """
    Parse GFA file once and extract coverage data for all specified contigs.
    Much more efficient than grepping for each contig individually.
    
    Args:
        gfa_file (str): Path to GFA file
        contig_ids (set): Set of contig IDs to extract data for
    
    Returns:
        dict: Dictionary mapping contig_id -> coverage_data tuple
    """
    print(f"Parsing GFA file {gfa_file} for {len(contig_ids)} contigs...")
    
    # Initialize data structures for each contig
    contig_data = {
        contig_id: {
            'positions': [],
            'dp1_values': [],
            'dp2_values': [],
            'dp3_values': [],
            'read_ids': [],
            'read_lens': [],
            'ol_lens': [],
            'snp_shares': [],
            'snp_diffs': [],
            'current_pos': 0
        }
        for contig_id in contig_ids
    }
    
    current_contig = None
    lines_processed = 0
    
    # Single pass through the GFA file
    with open(gfa_file, 'r') as f:
        for line in f:
            lines_processed += 1
            if lines_processed % 100000 == 0:
                print(f"  Processed {lines_processed:,} lines...")
            
            # Check for S lines (segment/contig headers) to identify which contig we're in
            if line.startswith('S'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    segment_id = parts[1]
                    # Check if this is one of our contigs
                    if segment_id in contig_ids:
                        current_contig = segment_id
                    else:
                        current_contig = None
                continue
            
            # Process alignment lines for the current contig
            if not line.startswith('a'):
                continue
            
            # If we're not tracking a contig, skip
            if current_contig is None:
                continue
            
            # Parse alignment line - use the same logic as plot.py
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            try:
                data = contig_data[current_contig]
                
                # Extract read ID from column 5 (index 4)
                read_id = parts[4] if len(parts) > 4 else "unknown_read"
                if len(read_id.split()) > 1:
                    read_id = read_id.split()[0] + ' ' + read_id.split()[-1]
                
                # Parse DP values and metrics from last column (same as plot.py)
                dp_part = parts[-1]
                metrics = {}
                for metric in dp_part.split(','):
                    if ':' in metric:
                        key, value = metric.split(':', 1)
                        try:
                            metrics[key] = int(value)
                        except ValueError:
                            metrics[key] = value
                
                # Extract values (same as plot.py)
                dp1 = metrics.get('DP1', 0)
                dp2 = metrics.get('DP2', 0)
                dp3 = metrics.get('DP3', 0)
                read_len = metrics.get('READ_LEN', 0)
                ol_len = metrics.get('OL_LEN_NEXT', 0)
                snp_share = metrics.get('SNP_SHARE_NEXT', 0)
                snp_diff = metrics.get('SNP_DIFF_NEXT', 0)
                
                # Update position (column 7, index 6)
                length = int(parts[6])
                data['current_pos'] += length
                
                # Store all data
                data['positions'].append(data['current_pos'])
                data['dp1_values'].append(dp1)
                data['dp2_values'].append(dp2)
                data['dp3_values'].append(dp3)
                data['read_ids'].append(read_id)
                data['read_lens'].append(read_len)
                data['ol_lens'].append(ol_len)
                data['snp_shares'].append(snp_share)
                data['snp_diffs'].append(snp_diff)
                
            except (IndexError, ValueError, KeyError) as e:
                print(f"Warning: Error parsing line for contig {current_contig}: {e}")
                continue
    
    print(f"Finished parsing GFA file ({lines_processed:,} lines)")
    
    # Convert to the format expected by create_interactive_plot
    result = {}
    for contig_id, data in contig_data.items():
        if data['positions']:  # Only include if we found data
            result[contig_id] = (
                data['positions'],
                data['dp1_values'],
                data['dp2_values'],
                data['dp3_values'],
                data['read_ids'],
                data['read_lens'],
                data['ol_lens'],
                data['snp_shares'],
                data['snp_diffs']
            )
            print(f"  Found {len(data['positions'])} alignment records for {contig_id}")
        else:
            print(f"  Warning: No coverage data found for {contig_id}")
            result[contig_id] = ([], [], [], [], [], [], [], [], [])
    
    return result


def calculate_contig_metrics(contig_info, coverage_data_dict, window_size=1000):
    """
    Calculate GC content and average coverage for each contig.
    
    Args:
        contig_info (dict): Dictionary with contig metadata (from extract_long_contigs)
        coverage_data_dict (dict): Dictionary with coverage data for each contig
        window_size (int): Window size for GC content calculation
    
    Returns:
        dict: Updated contig_info with 'gc_content', 'mean_coverage', and 'gc_data' added
    """
    for contig_id, info in contig_info.items():
        # Read contig sequence
        filepath = info['filepath']
        record = list(SeqIO.parse(filepath, "fasta"))[0]
        sequence = str(record.seq)
        
        # Calculate GC content across windows
        positions, gc_contents, gc_skew_contents = calculate_gc_content_windows(
            sequence, window_size
        )
        
        # Store GC data for later use in individual plots
        info['gc_data'] = (positions, gc_contents, gc_skew_contents)
        
        # Calculate mean GC content
        if gc_contents:
            info['gc_content'] = statistics.mean(gc_contents)
        else:
            info['gc_content'] = 0.0
        
        # Calculate mean coverage from DP1 values
        coverage_data = coverage_data_dict.get(contig_id)
        if coverage_data and coverage_data[1]:  # Check if dp1_values exist
            dp1_values = coverage_data[1]
            if dp1_values:
                info['mean_coverage'] = statistics.mean(dp1_values)
            else:
                info['mean_coverage'] = 0.0
        else:
            info['mean_coverage'] = 0.0
        
        print(f"  {contig_id}: GC={info['gc_content']:.2f}%, Coverage={info['mean_coverage']:.2f}x")
    
    return contig_info


def _generate_single_plot(args):
    """
    Helper function to generate a single plot (for parallel processing).
    
    Args:
        args: Tuple of (contig_id, info, coverage_data, plot_folder)
    
    Returns:
        Tuple of (contig_id, output_file)
    """
    contig_id, info, coverage_data, plot_folder = args
    
    # Use the full contig ID for the filename
    full_id = info['full_id']
    output_file = str(Path(plot_folder) / f"{full_id}.html")
    
    gc_data = info['gc_data']
    
    # Generate the plot
    create_interactive_plot(contig_id, gc_data, coverage_data, output_file)
    
    return contig_id, output_file


def generate_individual_plots(contig_info, coverage_data_dict, plot_folder, max_workers=5):
    """
    Generate individual interactive HTML plots for each contig in parallel.
    
    Args:
        contig_info (dict): Dictionary with contig metadata including gc_data
        coverage_data_dict (dict): Dictionary with coverage data for each contig
        plot_folder (str): Directory to store HTML plot files
        max_workers (int): Maximum number of parallel workers (default: 5)
    
    Returns:
        dict: Dictionary mapping contig_id -> plot_filepath
    """
    Path(plot_folder).mkdir(parents=True, exist_ok=True)
    
    # Limit workers to the number of contigs if fewer contigs than workers
    max_workers = min(max_workers, len(contig_info))
    
    print(f"  Generating plots using {max_workers} parallel workers...")
    
    # Prepare arguments for parallel processing
    plot_args = [
        (contig_id, info, coverage_data_dict[contig_id], plot_folder)
        for contig_id, info in contig_info.items()
    ]
    
    plot_paths = {}
    
    # Generate plots in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(_generate_single_plot, arg): arg[0] for arg in plot_args}
        
        # Collect results as they complete
        completed = 0
        total = len(futures)
        for future in as_completed(futures):
            contig_id = futures[future]
            try:
                result_contig_id, output_file = future.result()
                plot_paths[result_contig_id] = output_file
                completed += 1
                if completed % 5 == 0 or completed == total:
                    print(f"    Progress: {completed}/{total} plots completed")
            except Exception as e:
                print(f"    Error generating plot for {contig_id}: {e}")
    
    return plot_paths


def create_summary_report(contig_info, plot_paths, output_file, x_axis='gc'):
    """
    Create an interactive summary report showing all contigs.
    
    Args:
        contig_info (dict): Dictionary with contig metadata
        plot_paths (dict): Dictionary mapping contig_id -> plot file path
        output_file (str): Path for output HTML file
        x_axis (str): 'gc' for GC content or 'length' for contig length on x-axis
    """
    print(f"\nCreating summary report with x-axis={x_axis}...")

    # Some contigs may have failed plot generation (e.g. a plotting error) and
    # therefore have no entry in plot_paths. Exclude them from the report
    # instead of crashing.
    missing_plots = [cid for cid in contig_info if cid not in plot_paths]
    if missing_plots:
        print(f"  Warning: {len(missing_plots)} contig(s) have no plot and will be "
              f"excluded from the summary report: {', '.join(missing_plots)}")
        contig_info = {cid: info for cid, info in contig_info.items() if cid in plot_paths}

    # Separate circular and non-circular contigs
    circular_data = {'ids': [], 'lengths': [], 'gc': [], 'cov': [], 'links': [], 'sizes': [], 'hover': [], 'metadata': []}
    noncircular_data = {'ids': [], 'lengths': [], 'gc': [], 'cov': [], 'links': [], 'sizes': [], 'hover': [], 'metadata': []}
    
    all_lengths = []
    for contig_id, info in contig_info.items():
        all_lengths.append(info['length'])
    
    # Normalize sizes for better visualization
    min_size = 10
    max_size = 50
    if len(all_lengths) > 1:
        min_len = min(all_lengths)
        max_len = max(all_lengths)
    else:
        min_len = max_len = all_lengths[0] if all_lengths else 1
    
    for contig_id, info in contig_info.items():
        # Determine which dataset to add to
        target = circular_data if info['is_circular'] else noncircular_data
        
        target['ids'].append(contig_id)
        target['lengths'].append(info['length'])
        target['gc'].append(info['gc_content'])
        target['cov'].append(info['mean_coverage'])
        
        # Create relative path: plots/contig.html
        plot_filename = os.path.basename(plot_paths[contig_id])
        target['links'].append(f"plots/{plot_filename}")
        
        # Calculate size proportional to length
        if max_len > min_len:
            size = min_size + (max_size - min_size) * (info['length'] - min_len) / (max_len - min_len)
        else:
            size = max_size
        target['sizes'].append(size)
        
        # Format metadata for hover text (excluding 'len' since we show Length separately)
        metadata_lines = []
        for key, value in info['metadata'].items():
            if key.lower() != 'len':  # Skip the len field
                # Format the key-value pairs nicely
                formatted_key = key.replace('_', ' ').title()
                metadata_lines.append(f"{formatted_key}: {value}")
        
        # Create hover text with metadata
        hover_text = f"<b>{contig_id}</b><br>"
        hover_text += f"Length: {info['length']:,} bp<br>"
        hover_text += f"GC Content: {info['gc_content']:.2f}%<br>"
        hover_text += f"Mean Coverage: {info['mean_coverage']:.2f}x<br>"
        if metadata_lines:
            hover_text += "<br>" + "<br>".join(metadata_lines) + "<br>"
        hover_text += "<i>Click to view detailed plot</i>"
        
        target['hover'].append(hover_text)
    
    # Determine x-axis data and label
    if x_axis == 'length':
        x_label = 'Contig Length (bp)'
        x_circular = circular_data['lengths']
        x_noncircular = noncircular_data['lengths']
    else:  # x_axis == 'gc'
        x_label = 'Mean GC Content (%)'
        x_circular = circular_data['gc']
        x_noncircular = noncircular_data['gc']
    
    # Create the scatter plot with separate traces for circular and non-circular
    fig = go.Figure()
    
    # Add non-circular contigs (squares) - with both x-axis data stored
    if noncircular_data['ids']:
        fig.add_trace(go.Scatter(
            x=x_noncircular,
            y=noncircular_data['cov'],
            mode='markers',
            marker=dict(
                size=noncircular_data['sizes'],
                color=noncircular_data['lengths'],
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(
                    title='Length (bp)',
                    thickness=15,
                    len=0.7
                ),
                symbol='square',
                line=dict(width=1, color='black')
            ),
            text=noncircular_data['ids'],
            hovertext=noncircular_data['hover'],
            hoverinfo='text',
            customdata=[[link, gc, length] for link, gc, length in zip(
                noncircular_data['links'], 
                noncircular_data['gc'], 
                noncircular_data['lengths']
            )],
            name='Linear Contigs',
            visible=True
        ))
    
    # Add circular contigs (circles) - with both x-axis data stored
    if circular_data['ids']:
        fig.add_trace(go.Scatter(
            x=x_circular,
            y=circular_data['cov'],
            mode='markers',
            marker=dict(
                size=circular_data['sizes'],
                color=circular_data['lengths'],
                colorscale='Viridis',
                showscale=False,  # Only show one colorbar
                symbol='circle',
                line=dict(width=1, color='black')
            ),
            text=circular_data['ids'],
            hovertext=circular_data['hover'],
            hoverinfo='text',
            customdata=[[link, gc, length] for link, gc, length in zip(
                circular_data['links'], 
                circular_data['gc'], 
                circular_data['lengths']
            )],
            name='Circular Contigs',
            visible=True
        ))
    
    # Update layout with interactive controls
    fig.update_layout(
        title=dict(
            text=f'Contig Summary Report<br><sub>{len(contig_info)} contigs analyzed</sub>',
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=x_label,
            side='bottom'
        ),
        yaxis_title='Mean Coverage (x)',
        width=1000,
        height=750,
        hovermode='closest',
        template='plotly_white',
        font=dict(size=12),
        updatemenus=[
            # Y-axis scale toggle
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{"yaxis.type": "linear"}],
                        label="Linear Y",
                        method="relayout"
                    ),
                    dict(
                        args=[{"yaxis.type": "log"}],
                        label="Log Y",
                        method="relayout"
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top"
            ),
            # X-axis toggle (moved below Y-axis buttons)
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        label="X: GC Content",
                        method="update",
                        args=[
                            {"x": [noncircular_data['gc'] if noncircular_data['gc'] else [], 
                                   circular_data['gc'] if circular_data['gc'] else []]},
                            {"xaxis.title.text": "Mean GC Content (%)"}
                        ]
                    ),
                    dict(
                        label="X: Length",
                        method="update",
                        args=[
                            {"x": [noncircular_data['lengths'] if noncircular_data['lengths'] else [], 
                                   circular_data['lengths'] if circular_data['lengths'] else []]},
                            {"xaxis.title.text": "Contig Length (bp)"}
                        ]
                    )
                ],
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.06,
                yanchor="top"
            ),
        ]
    )
    
    # Add JavaScript for click handling to open individual plots
    fig.update_traces(
        marker=dict(opacity=0.8),
        selector=dict(mode='markers')
    )
    
    # Save the plot with custom JavaScript for click handling
    html_string = plot(fig, output_type='div', include_plotlyjs='cdn', config={'displayModeBar': True})
    
    # Prepare contig list for dropdown (sorted by name)
    import json
    all_contigs = []
    for contig_id, info in contig_info.items():
        all_contigs.append({
            'id': contig_id,
            'full_id': info['full_id'],
            'link': f"plots/{os.path.basename(plot_paths[contig_id])}",
            'length': info['length'],
            'is_circular': info['is_circular']
        })
    all_contigs.sort(key=lambda x: x['id'])
    contig_data_json = json.dumps(all_contigs)
    
    # Add custom JavaScript for click handling
    custom_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Contig Summary Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            text-align: center;
            margin-bottom: 20px;
        }}
        .info {{
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .controls {{
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .control-group {{
            margin-bottom: 15px;
        }}
        .control-group label {{
            display: inline-block;
            width: 150px;
            font-weight: bold;
        }}
        .control-group input[type="number"] {{
            width: 150px;
            padding: 5px;
            border: 1px solid #ccc;
            border-radius: 3px;
        }}
        .control-group input[type="checkbox"] {{
            margin-left: 10px;
        }}
        .control-group select {{
            width: 300px;
            padding: 5px;
            border: 1px solid #ccc;
            border-radius: 3px;
        }}
        .control-group button {{
            padding: 8px 15px;
            background-color: #007bff;
            color: white;
            border: none;
            border-radius: 3px;
            cursor: pointer;
            margin-left: 10px;
        }}
        .control-group button:hover {{
            background-color: #0056b3;
        }}
        .plot-container {{
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
    </style>
</head>
<body>
    <div class="info">
        <h3>Instructions</h3>
        <ul>
            <li><b>Hover</b> over a point to see contig details and metadata</li>
            <li><b>Click</b> on a point to open the detailed plot for that contig</li>
            <li><b>Dropdown menu</b>: Search for and select a contig by name to open its plot</li>
            <li><b>Circles</b> represent circular contigs</li>
            <li><b>Squares</b> represent linear contigs</li>
            <li>Point size is proportional to contig length</li>
            <li>Point color represents contig length (see colorbar)</li>
            <li>Use the <b>Linear Y / Log Y</b> and <b>X: GC / X: Length</b> buttons above the plot to toggle axes</li>
        </ul>
    </div>
    <div class="controls">
        <h3>Filters</h3>
        <div class="control-group">
            <label for="minLength">Minimum Length (bp):</label>
            <input type="number" id="minLength" value="0" min="0" step="10000">
            <button onclick="applyFilters()">Apply Filter</button>
            <button onclick="resetFilters()">Reset</button>
        </div>
        <div class="control-group">
            <label for="circularOnly">Show only circular:</label>
            <input type="checkbox" id="circularOnly" onchange="applyFilters()">
        </div>
        <div class="control-group">
            <label for="contigSelect">Jump to contig:</label>
            <select id="contigSelect" onchange="openContigPlot()">
                <option value="">-- Select a contig --</option>
            </select>
        </div>
    </div>
    <div class="plot-container">
        {html_string}
    </div>
    <script>
        // Contig data
        var contigData = {contig_data_json};
        
        // Populate contig dropdown
        var select = document.getElementById('contigSelect');
        contigData.forEach(function(contig) {{
            var option = document.createElement('option');
            option.value = contig.link;
            var label = contig.id;
            if (contig.is_circular) {{
                label += ' (circular)';
            }}
            label += ' - ' + contig.length.toLocaleString() + ' bp';
            option.textContent = label;
            select.appendChild(option);
        }});
        
        // Function to open selected contig plot
        function openContigPlot() {{
            var select = document.getElementById('contigSelect');
            var link = select.value;
            if (link) {{
                window.open(link, '_blank');
                select.value = '';  // Reset selection
            }}
        }}
        
        // Store original marker sizes for reset
        var originalSizes = [];
        
        // Add click handler to scatter plot
        var myPlot = document.getElementsByClassName('plotly-graph-div')[0];
        myPlot.on('plotly_click', function(data){{
            var point = data.points[0];
            if (point.customdata && point.customdata[0]) {{
                var plotPath = point.customdata[0];
                window.open(plotPath, '_blank');
            }}
        }});
        
        // Store original sizes after plot loads
        setTimeout(function() {{
            var data = myPlot.data;
            for (var i = 0; i < data.length; i++) {{
                if (Array.isArray(data[i].marker.size)) {{
                    originalSizes[i] = data[i].marker.size.slice();
                }} else {{
                    originalSizes[i] = data[i].marker.size;
                }}
            }}
        }}, 100);
        
        // Filter functions
        function applyFilters() {{
            var minLength = parseInt(document.getElementById('minLength').value) || 0;
            var circularOnly = document.getElementById('circularOnly').checked;
            
            // Get current data from the plot
            var data = myPlot.data;
            
            for (var i = 0; i < data.length; i++) {{
                var trace = data[i];
                
                // Determine which trace this is based on name
                var isCircular = trace.name === 'Circular Contigs';
                
                // If circular only is checked, hide linear contigs entirely
                if (circularOnly && !isCircular) {{
                    Plotly.restyle(myPlot, {{'visible': false}}, [i]);
                    continue;
                }} else {{
                    Plotly.restyle(myPlot, {{'visible': true}}, [i]);
                }}
                
                // Filter by length using customdata[2] which contains length
                if (trace.customdata && originalSizes[i]) {{
                    var newSizes;
                    if (Array.isArray(originalSizes[i])) {{
                        newSizes = originalSizes[i].map(function(size, idx) {{
                            var length = trace.customdata[idx][2];  // d[2] is length
                            return length >= minLength ? size : 0;  // Set size to 0 to hide
                        }});
                    }} else {{
                        newSizes = originalSizes[i];
                    }}
                    
                    Plotly.restyle(myPlot, {{'marker.size': [newSizes]}}, [i]);
                }}
            }}
        }}
        
        function resetFilters() {{
            document.getElementById('minLength').value = 0;
            document.getElementById('circularOnly').checked = false;
            
            // Reset all traces to visible with original sizes
            var data = myPlot.data;
            for (var i = 0; i < data.length; i++) {{
                if (originalSizes[i]) {{
                    Plotly.restyle(myPlot, {{
                        'marker.size': [originalSizes[i]],
                        'visible': true
                    }}, [i]);
                }}
            }}
        }}
    </script>
</body>
</html>
"""
    
    # Write the HTML file
    with open(output_file, 'w') as f:
        f.write(custom_html)
    
    print(f"Summary report saved as {output_file}")


def main(args):
    """
    Main function for report generation.
    
    Args:
        args: Parsed command-line arguments
    """
    print(f"\n{'='*60}")
    print(f"MyloTools Report Generator")
    print(f"{'='*60}\n")
    
    # Step 1: Extract long contigs
    print(f"Step 1: Extracting contigs > {args.min_length:,} bp...")
    contig_folder = os.path.join(args.output, 'contigs')
    contig_info = extract_long_contigs(args.fasta, args.min_length, contig_folder)
    print(f"  Extracted {len(contig_info)} contigs\n")
    
    if not contig_info:
        print("No contigs found meeting the length threshold. Exiting.")
        return
    
    # Step 2: Parse GFA file efficiently
    print(f"Step 2: Parsing GFA file for coverage data...")
    contig_ids = set(contig_info.keys())
    coverage_data_dict = parse_gfa_for_all_contigs(args.gfa, contig_ids)
    print()
    
    # Step 3: Calculate metrics (GC content and mean coverage)
    print(f"Step 3: Calculating contig metrics...")
    contig_info = calculate_contig_metrics(
        contig_info, coverage_data_dict, args.window_size
    )
    print()
    
    # Step 4: Generate individual plots
    print(f"Step 4: Generating individual contig plots...")
    plot_folder = os.path.join(args.output, 'plots')
    plot_paths = generate_individual_plots(
        contig_info, 
        coverage_data_dict, 
        plot_folder, 
        max_workers=args.workers
    )
    print(f"  Generated {len(plot_paths)} plots\n")
    
    # Step 5: Create summary report
    print(f"Step 5: Creating summary report...")
    summary_file = os.path.join(args.output, 'contig_summary_report.html')
    create_summary_report(contig_info, plot_paths, summary_file, args.x_axis)
    
    print(f"\n{'='*60}")
    print(f"Report generation complete!")
    print(f"{'='*60}")
    print(f"\nOutput files:")
    print(f"  - Extracted contigs: {contig_folder}/")
    print(f"  - Individual plots: {plot_folder}/")
    print(f"  - Summary report: {summary_file}")
    print(f"\nOpen {summary_file} in a web browser to explore your assembly.\n")
