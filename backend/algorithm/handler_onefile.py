#!/usr/bin/env python3
"""
CSV to JSON converter for spline lab tables.

Usage:
    python handler_onefile.py [input.csv]
    If no file is given, reads from standard input.

The script automatically detects the table type
(coefficients, comparison, or convergence) and outputs
a corresponding JSON structure.
"""

import sys
import csv
import json
from typing import List

# Recognisable headers (order matters)
TABLE1_HEADER = ['i', 'x_{i-1}', 'x_i', 'a', 'b', 'c', 'd']
TABLE2_HEADER = [
    'j', 'x', 'F(x)', 'S(x)', 'diff_F',
    "F'(x)", "S'(x)", "diff_F'",
    "F''(x)", "S''(x)", "diff_F''"
]
CONV_HEADER = ['n', 'max_err_F', "max_err_F'", "max_err_F''"]


def detect_table(headers: List[str]) -> str:
    """Return table type based on the headers."""
    h = [s.strip() for s in headers]
    if h == TABLE1_HEADER:
        return 'table1'
    if h == TABLE2_HEADER:
        return 'table2'
    if h == CONV_HEADER:
        return 'convergence'
    return 'unknown'


def parse_row(row: List[str], table_type: str) -> dict:
    """Convert a CSV row into a dictionary depending on the table type."""
    try:
        if table_type == 'table1':
            return {
                "i": int(row[0]),
                "x_{i-1}": float(row[1]),
                "x_i": float(row[2]),
                "a": float(row[3]),
                "b": float(row[4]),
                "c": float(row[5]),
                "d": float(row[6])
            }
        if table_type == 'table2':
            return {
                "j": int(row[0]),
                "x": float(row[1]),
                "F": float(row[2]),
                "S": float(row[3]),
                "diff_F": float(row[4]),
                "F'": float(row[5]),
                "S'": float(row[6]),
                "diff_F'": float(row[7]),
                "F''": float(row[8]),
                "S''": float(row[9]),
                "diff_F''": float(row[10])
            }
        if table_type == 'convergence':
            return {
                "n": int(row[0]),
                "max_err_F": float(row[1]),
                "max_err_F'": float(row[2]),
                "max_err_F''": float(row[3])
            }
        # Fallback for unknown formats
        return {f"col{i}": val for i, val in enumerate(row)}
    except (ValueError, IndexError) as e:
        return {"error": f"malformed row: {e}", "raw": row}


def main() -> None:
    # Choose input source
    if len(sys.argv) > 1:
        infile = open(sys.argv[1], 'r', newline='')
    else:
        infile = sys.stdin

    reader = csv.reader(infile)
    try:
        headers = next(reader)
    except StopIteration:
        print(json.dumps({"error": "empty input"}))
        return

    table_type = detect_table(headers)
    data = []
    for row in reader:
        if not row:               # skip blank lines
            continue
        data.append(parse_row(row, table_type))

    result = {
        "type": table_type,
        "data": data
    }
    # Print JSON to stdout
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()