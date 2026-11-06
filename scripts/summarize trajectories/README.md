# Summarize identified trajectories

The scripts and the README.md in this folder for summarizing the trajectories identified by POTTR were generated using Claude Code. However, we, the authors, reviewed the code carefully. The code in this directory was only used to create the overview visualizations of our identified trajectories. We included the visualizations the supplemental results.

Draws every trajectory reported by `get_output_sizes.py` (one `run_POTTR.py` output folder per k) with its support, supporting trees, p-values and resolution evidence, and merges everything into one pdf.

Each trajectory appears only once. If the same trajectory was reported for several k, the occurrence with the highest support is shown (ties: more supporting trees, then one with p-values, then the smallest k). Two trajectories count as the same if their alterations, clusters and transitively closed relations are identical. Trajectories are sorted by support, then by size (both descending).

The drawings are placed on A4 pages in two columns (left to right, top to bottom). A drawing that would be scaled below 60 % in a column gets a row over the full page width.

The values shown above each drawing are read from `out_k<k>/significance_output.txt`:
- `pval_ind`: MASTRO statistical test (independent model)
- `pval_POTTR`: POTTR statistical test
- `resolution_evidence`: evidence score for cluster resolution (1 = all supporting patients fit the trajectory without resolving clusters, smaller = more orders rely on resolved clusters)

`run_POTTR.py` only runs the significance test for trajectories with fewer than 13 nodes (root included). For larger trajectories, and for older `significance_output.txt` files without a column, the summary shows `n/a`.

## Usage

```bash
cd Code
python "summarize trajectories/summarize_trajectories.py" \
    --results <output dir of get_output_sizes.py, e.g. .../2026-10-01_12-00-00/> \
    --output <summary output dir>
```

Options:
- `--max-trajectories-per-k N`: only read the first N trajectories per k (useful with large solution pools).
- `--labels-per-line N`: cluster members per line in a node label (default 4).
- `--name`: base name of the pdf and csv (default `tracerx_trajectory_summary`).

Output in `--output`:
- `<name>.pdf`: the trajectories on A4 pages, with a bookmark per trajectory. Drawings wider or higher than their column or page are scaled down.
- `<name>.csv`: one row per distinct trajectory: number, size, support, supporting trees, `pval_ind`, `pval_POTTR`, `resolution_evidence`, note, all k it was found for, pdf page, relations.
- `drawings/`: the single `.dot` and `.pdf` drawings.

## Files
- `parse_results.py`: finds the `out_k*` folders, reads `converted_graphs.txt` and `significance_output.txt` (reuses `create_graphs.get_graph_from_line`) and merges trajectories found for several k.
- `draw_graphs.py`: Graphviz drawings. Clusters (`-?-`) are merged into one grey box, and edges are drawn as a transitive reduction with root `0`.
- `create_overview.py`: csv with one row per trajectory.
- `summarize_trajectories.py`: command line script that puts everything together and places the drawings on A4 pages.

## Requirements
The POTTR environment (networkx), Graphviz (`dot` on the PATH) and `pypdf` (`pip install pypdf`).
