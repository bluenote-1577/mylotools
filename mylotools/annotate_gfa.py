#!/usr/bin/env python3
"""
Annotate a myloasm final_contig_graph.gfa with real sequence bases from the
final assembly.

The GFA emitted by myloasm is written before polishing and before contig
filtering, so `S` line lengths can differ from the final assembly and some
GFA segments have no corresponding contig in the final assembly at all. This
fills in each segment's sequence field with the matching contig's bases from
assembly_primary.fa, and marks segments that were filtered out of the final
assembly with `*`.
"""
from Bio import SeqIO


def load_assembly_sequences(fasta_file):
    """Map base contig ID (e.g. 'u273ctg') -> sequence string."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        base_id = record.id.split('_')[0]
        sequences[base_id] = str(record.seq)
    return sequences


def annotate_gfa(gfa_file, fasta_file, output_file, keep_read_info=False):
    """
    Write an annotated copy of gfa_file to output_file.

    Returns:
        dict with counts: annotated, filtered, read_lines_dropped
    """
    sequences = load_assembly_sequences(fasta_file)

    counts = {'annotated': 0, 'filtered': 0, 'read_lines_dropped': 0}

    with open(gfa_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('a') and not keep_read_info:
                counts['read_lines_dropped'] += 1
                continue

            if not line.startswith('S'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            segment_id = parts[1]
            seq = sequences.get(segment_id)

            if seq is None:
                parts[2] = '*'
                parts.append('FILTERED:Z:FAIL')
                counts['filtered'] += 1
            else:
                parts[2] = seq
                for i, field in enumerate(parts):
                    if field.startswith('LN:i:'):
                        gfa_length = field[len('LN:i:'):]
                        parts[i] = f'LN:i:{len(seq)}'
                        parts.append(f'LN_GFA:i:{gfa_length}')
                        break
                parts.append('FILTERED:Z:GOOD')
                counts['annotated'] += 1

            fout.write('\t'.join(parts) + '\n')

    return counts


def main(args):
    print(f"\n{'='*60}")
    print(f"Mylotools GFA Annotator")
    print(f"{'='*60}\n")

    print(f"Loading final assembly sequences from {args.fasta}...")
    print(f"Annotating {args.gfa} -> {args.output}...")

    counts = annotate_gfa(
        args.gfa, args.fasta, args.output,
        keep_read_info=args.keep_read_info
    )

    print(f"\n{'='*60}")
    print(f"Annotation complete!")
    print(f"{'='*60}")
    print(f"  Segments annotated with bases: {counts['annotated']}")
    print(f"  Segments left as '*' (not in final assembly): {counts['filtered']}")
    if not args.keep_read_info:
        print(f"  Read-alignment ('a') lines dropped: {counts['read_lines_dropped']}")
    print(f"\nWrote {args.output}\n")
