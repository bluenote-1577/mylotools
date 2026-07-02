# Changelog

All notable changes to mylotools are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [2.1.0]

### Added
- `mylotools annotate-gfa`: annotate a myloasm `final_contig_graph.gfa` with real sequence bases from the final assembly fasta.
  - The GFA myloasm emits is written before polishing and before contig filtering, so `S` line lengths can differ from the final assembly and some GFA segments don't exist in the final assembly at all.
  - Segments matching a contig in the final assembly get their `*` sequence field filled in with actual bases, `LN:i:` updated to the polished length, the original pre-polish length preserved as `LN_GFA:i:`, and tagged `FILTERED:Z:GOOD`.
  - Segments with no match in the final assembly (filtered out post-assembly) are left as `*` and tagged `FILTERED:Z:FAIL`.
  - Per-read alignment (`a`) lines are dropped by default, since their coordinates refer to the pre-polish sequence; pass `--keep-read-info` to retain them.
  - `L` (link) lines are left untouched.

## [2.0.0]

### Added
- `mylotools report`: generate a comprehensive HTML report with per-contig plots (GC content, GC skew, coverage, read overlaps) and a summary scatter plot for all long/circular contigs.
- `mylotools sanitize-headers`: sanitize FASTA headers by replacing underscores with spaces (except the first), with automatic backup.
