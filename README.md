# mylotools - utility scripts for visualizing/manipulating myloasm's outputs.

These are utility scripts for manipulating the output for [myloasm](https://github.com/bluenote-1577/myloasm). 

The documentation is hosted at [the myloasm manual](https://myloasm-docs.github.io/). 

See [CHANGELOG.md](CHANGELOG.md) for release history.

## Commands

- `report` - generate a comprehensive HTML report with per-contig plots and a summary scatter plot for all long/circular contigs.
- `plot` - generate a plot of various statistics (GC content, coverage, read overlaps) for one contig.
- `strain-viz` - visualize overlaps between and within two or more similar contigs.
- `extract-contigs` - extract contigs longer than a given length into their own fasta files.
- `sanitize-headers` - replace underscores with spaces in FASTA headers (keeps a backup).
- `annotate-gfa` - fill in real sequence bases in `final_contig_graph.gfa` from the final assembly.

### `annotate-gfa`

myloasm's `final_contig_graph.gfa` is written **before polishing and before contig filtering**, so it disagrees
with the final assembly fasta (e.g. `assembly_primary.fa`) in two ways: segment lengths differ from the polished
contigs, and some GFA segments don't survive into the final assembly at all. This command reconciles the two by
copying real sequence bases from the final assembly into the GFA's `S` lines.

```
mylotools annotate-gfa --gfa final_contig_graph.gfa --fasta assembly_primary.fa --output annotated.gfa
```

- If no match is found (the contig was filtered out after assembly), the sequence field is left as `*` and the segment is tagged `FILTERED:Z:FAIL`.

- `L` (link) lines are left untouched, so the graph structure stays valid regardless of which segments were filtered.

- By default, per-read alignment (`a`) lines are dropped from the output, since their coordinates refer to the pre-polish sequence and no longer line up once bases are swapped in. Pass `--keep-read-info` to retain them as-is.

