#!/usr/bin/env python3
"""Convert mean_variance_change_summary.txt into a markdown comparison table row.

Usage:
    python3 make_table.py <summary.txt> [--header|--no-header]

If --header is passed (or no flag at all), the header + separator lines are
printed above the data row. If --no-header is passed, only the data row is
printed (useful for stitching multiple experiments into one table).
"""

import sys
import re

# metric name -> short label used in the table header
LABELS = {
    "kmer_count": "kmer",
    "color_class_count": "CT",
    "overall_index_size": "OSize",
    "color_class_table_size": "CTSZ",
    "mst_color_class_table_size": "MSTCT",
}

# Order metrics should appear in the output table
ORDER = ["kmer_count", "color_class_count", "overall_index_size",
         "color_class_table_size", "mst_color_class_table_size"]


def parse_file(path, mean:bool = True):
    offset = 0 if mean else 2
    data = {}
    with open(path) as f:
        next(f)  # skip header line
        for line in f:
            parts = line.split()
            if not parts:
                continue
            metric = parts[0]
            if metric not in LABELS:
                continue  # skip file_count, cqf_size, or anything unrecognized
            mean_change = float(parts[1+offset])
            variance_change = float(parts[2+offset])
            data[metric] = (mean_change, variance_change)
    return data


def build_headers(exp_mode:bool = 1):
    headers = ["# Experiments"] if exp_mode else ["# Workers"]
    for metric in ORDER:
        label = LABELS[metric]
        headers.append(f"{label}<br>RR")
        headers.append(f"{label}<br>WR")
    return headers


def build_row(data, n_experiments):
    row = [str(n_experiments)]
    for metric in ORDER:
        rr, wr = data.get(metric, (0.0, 0.0))
        row.append(f"{rr:.3f}%")
        row.append(f"{wr:.3f}%")
    return row


def fmt_row(values, widths):
    return "| " + " | ".join(v.ljust(w) for v, w in zip(values, widths)) + " |"


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 make_table.py <summary.txt> [--header|--no-header]")
        sys.exit(1)

    path = sys.argv[1]

    # parse the header flag (defaults to printing the header)
    print_header = False
    if len(sys.argv) >= 3:
        flag = sys.argv[2].lower()
        if flag in ("--no-header", "no-header", "false", "0"):
            print_header = False
        elif flag in ("--header", "header", "true", "1"):
            print_header = True
        else:
            print(f"Unrecognized flag: {sys.argv[2]!r}. Use --header or --no-header.")
            sys.exit(1)
    exp_mode = True
    if len(sys.argv) >= 4:
        flag = sys.argv[3].lower()
        if flag in ("--experiments", "-e", "false", "0"):
            exp_mode = True
        elif flag in ("--workers", "-w", "true", "1"):
            exp_mode = False
        else:
            print(f"Unrecognized flag: {sys.argv[3]!r}")
            sys.exit(1)
    mean = True
    if len(sys.argv) >= 5:
        flag = sys.argv[4].lower()
        if flag in ("--mean", "-m", "true", "1"):
            mean = True
        elif flag in ("--variance", "-v", "false", "0"):
            mean = False
        else:
            print(f"Unrecognized flag: {sys.argv[4]!r}")
            sys.exit(1)

    # infer n_experiments from the path (e.g. 500_sample_ids...)
    m = re.search(r"(\d+)_sample", path) if exp_mode else re.search(r"(\d+)_clusters", path)
    n_experiments = m.group(1) if m else "?"

    data = parse_file(path, mean)
    headers = build_headers(exp_mode)
    row = build_row(data, n_experiments)

    cols = list(zip(headers, row))
    widths = [max(len(h), len(v)) for h, v in cols]

    lines = []
    if print_header:
        lines.append(fmt_row(headers, widths))
        lines.append("| " + " | ".join("-" * w for w in widths) + " |")
    lines.append(fmt_row(row, widths))

    print("\n".join(lines))


if __name__ == "__main__":
    main()
