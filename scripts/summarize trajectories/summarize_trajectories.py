"""
Summarize all trajectories found by get_output_sizes.py (e.g. for the TRACERx data) in a single pdf.

Trajectories are collected from out_k<k>/converted_graphs.txt of every k. A trajectory that was reported for
several k is shown only once, with its highest support. Each trajectory is drawn together with its support, its
supporting trees, the p-value of the MASTRO statistical test (pval_ind, independent model), the p-value of the
POTTR statistical test (pval_POTTR) and the evidence score for cluster resolution (resolution_evidence), all read
from out_k<k>/significance_output.txt. The drawings are placed in two columns on A4 pages; drawings that would
become too small in a column span the full page width.

The pdf has a bookmark per trajectory. A csv with one row per trajectory is written next to the pdf.

Example:
    python "summarize trajectories/summarize_trajectories.py" \
        --results ../Output/output_sizes2026-10-01_12-00-00/ \
        --output ../Output/summary/

Requires Graphviz ('dot' on the PATH) and pypdf (pip install pypdf) in addition to the POTTR environment.
"""
import argparse
import os
import sys
import textwrap

from parse_results import load_results, unique_trajectories
from draw_graphs import to_dot, message_dot, render_pdf
from create_overview import write_csv, format_pval

# A4 in pt
A4_WIDTH = 595.28
A4_HEIGHT = 841.89
PAGE_MARGIN = 28        # ~1 cm
GAP = 10                # vertical space between two rows of trajectories
COLUMN_GAP = 14         # horizontal space between the left and right column
MIN_COLUMN_SCALE = 0.6  # drawings that would be scaled down more than this in a column span both columns
TITLE_FONTSIZE = 10
TITLE_WIDTH = 44        # characters per title line before wrapping


def get_parser():
    parser = argparse.ArgumentParser(description='Draw all trajectories reported by get_output_sizes.py '
                                                 'and their p-values in a single pdf.')
    parser.add_argument('--results', '-r', required=True, dest='results', type=str,
                        help='Output directory of get_output_sizes.py (contains out_k2, out_k3, ...) '
                             'or a single run_POTTR.py output directory')
    parser.add_argument('--output', '-o', required=True, dest='output', type=str,
                        help='Directory to store the summary pdf, csv and single drawings')
    parser.add_argument('--name', '-n', dest='name', type=str, default='trajectory_summary',
                        help='Base name of the summary pdf and csv (default trajectory_summary)')
    parser.add_argument('--max-trajectories-per-k', dest='max_per_k', type=int, default=None,
                        help='Only read the first N trajectories per k (default all)')
    parser.add_argument('--labels-per-line', dest='per_line', type=int, default=4,
                        help='Number of cluster members per line in a node label (default 4)')
    parser.add_argument('--dot', dest='dot', type=str, default='dot',
                        help='Graphviz dot executable (default dot)')
    return parser


def wrap(text: str):
    # do not break tree names like CRUK0063-0 at the hyphen
    return textwrap.wrap(text, width=TITLE_WIDTH, subsequent_indent='    ', break_on_hyphens=False) or [text]


def trajectory_title(t, n_trajectories):
    lines = [f'Trajectory {t.number} of {n_trajectories}  |  {t.size} alterations',
             f'support: {t.support} patients ({len(t.supporting_graphs)} trees)',
             f'MASTRO test, pval_ind: {format_pval(t.pval_ind)}',
             f'POTTR test, pval_POTTR: {format_pval(t.pval_pottr)}',
             f'resolution score: {format_pval(t.evidence)}']
    if t.significance_note:
        lines += wrap(f'note: {t.significance_note}')
    lines += wrap('supporting trees: ' + (', '.join(t.supporting_graphs) or '-'))
    return lines


def import_pypdf():
    try:
        import pypdf
    except ImportError:
        sys.exit('pypdf is required to merge the drawings into one pdf: pip install pypdf')
    return pypdf


class Summary:
    def __init__(self, args):
        self.args = args
        self.drawings_dir = os.path.join(args.output, 'drawings')
        self.entries = []      # (pdf file, outline title, trajectory number)

    def build(self, unique):
        for t in unique:
            title = trajectory_title(t, len(unique))
            if t.graph is None:
                dot_source = message_dot(title + ['could not build the trajectory graph'], TITLE_FONTSIZE)
            else:
                dot_source = to_dot(t.graph, title, per_line=self.args.per_line, title_fontsize=TITLE_FONTSIZE)
            pdf = render_pdf(dot_source, os.path.join(self.drawings_dir, f'trajectory_{t.number}.pdf'),
                             self.args.dot)
            outline = f'Trajectory {t.number} (size {t.size}, support {t.support})'
            self.entries.append((pdf, outline, t.number))

    def run(self):
        args = self.args
        os.makedirs(self.drawings_dir, exist_ok=True)

        results = load_results(args.results, args.max_per_k)
        if not results:
            sys.exit(f'No run_POTTR.py output folders (out_k<k>/converted_graphs.txt) found in {args.results}')
        for k, trajectories in sorted(results.items()):
            print(f'k = {k}: {len(trajectories)} trajectories')

        unique = unique_trajectories(results)
        print(f'{len(unique)} distinct trajectories')

        self.build(unique)
        placements = layout(self.entries)

        # page numbers are 1-based as shown by pdf viewers
        page_of = {number: 1 + p.page for (_, _, number), p in zip(self.entries, placements)}
        csv_file = write_csv(unique, page_of, os.path.join(args.output, args.name + '.csv'))

        out_pdf = os.path.join(args.output, args.name + '.pdf')
        merge(self.entries, placements, out_pdf)
        print(f'Summary pdf written to {out_pdf}')
        print(f'Trajectory table written to {csv_file}')


class Placement:
    def __init__(self, page, x, y, scale, top):
        self.page = page        # 0-based index among the drawing pages
        self.x = x              # lower left corner of the drawing on the A4 page (pt)
        self.y = y
        self.scale = scale
        self.top = top          # upper edge of the drawing, target of its bookmark


def layout(entries):
    """
    Place the drawings in rows of two columns (left to right, top to bottom) on A4 pages and start a new page
    when the next row does not fit. A drawing that would have to be scaled below MIN_COLUMN_SCALE to fit a column
    gets a row of its own over the full page width. Drawings are never enlarged.
    """
    pypdf = import_pypdf()
    full_width = A4_WIDTH - 2 * PAGE_MARGIN
    column_width = (full_width - COLUMN_GAP) / 2
    max_height = A4_HEIGHT - 2 * PAGE_MARGIN

    # (scaled width, scaled height, scale, spans both columns)
    items = []
    for pdf, _, _ in entries:
        box = pypdf.PdfReader(pdf).pages[0].mediabox
        width, height = float(box.width), float(box.height)
        scale = min(1.0, column_width / width, max_height / height)
        full_scale = min(1.0, full_width / width, max_height / height)
        wide = scale < MIN_COLUMN_SCALE and full_scale > scale
        if wide:
            scale = full_scale
        items.append((width * scale, height * scale, scale, wide))

    rows = []
    for i, (_, _, _, wide) in enumerate(items):
        if wide or not rows or len(rows[-1]) == 2 or items[rows[-1][0]][3]:
            rows.append([i])
        else:
            rows[-1].append(i)

    placements = [None] * len(items)
    top = A4_HEIGHT - PAGE_MARGIN
    page, y = 0, top
    for row in rows:
        row_height = max(items[i][1] for i in row)
        if y < top and y - row_height < PAGE_MARGIN:
            page, y = page + 1, top
        for column, i in enumerate(row):
            _, height, scale, _ = items[i]
            x = PAGE_MARGIN + column * (column_width + COLUMN_GAP)
            placements[i] = Placement(page, x, y - height, scale, y)
        y -= row_height + GAP
    return placements


def merge(entries, placements, out_pdf):
    pypdf = import_pypdf()
    from pypdf.generic import Fit

    writer = pypdf.PdfWriter()
    n_pages = max(p.page for p in placements) + 1 if placements else 0
    pages = [writer.add_blank_page(A4_WIDTH, A4_HEIGHT) for _ in range(n_pages)]

    for (pdf, outline, _), p in zip(entries, placements):
        drawing = pypdf.PdfReader(pdf).pages[0]
        box = drawing.mediabox
        transformation = (pypdf.Transformation()
                          .translate(-float(box.left), -float(box.bottom))
                          .scale(p.scale, p.scale)
                          .translate(p.x, p.y))
        pages[p.page].merge_transformed_page(drawing, transformation)
        writer.add_outline_item(outline, p.page, fit=Fit.xyz(left=p.x, top=p.top))

    with open(out_pdf, 'wb') as f:
        writer.write(f)


if __name__ == '__main__':
    Summary(get_parser().parse_args()).run()
