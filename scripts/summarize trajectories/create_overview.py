"""
Csv with one row per unique trajectory.
"""
import csv


def format_pval(p):
    return 'n/a' if p is None else f'{p:.3g}'


def write_csv(unique: list, page_of: dict, out_csv: str):
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['trajectory', 'size', 'support', 'supporting_trees', 'pval_ind', 'pval_POTTR',
                         'resolution_evidence', 'note', 'found_for_k', 'pdf_page', 'relations'])
        for t in unique:
            writer.writerow([t.number, t.size, t.support, ' '.join(t.supporting_graphs),
                             '' if t.pval_ind is None else t.pval_ind,
                             '' if t.pval_pottr is None else t.pval_pottr,
                             '' if t.evidence is None else t.evidence,
                             t.significance_note, ' '.join(map(str, t.found_for_k)),
                             page_of.get(t.number, ''), ' '.join(t.relations)])
    return out_csv
