import argparse
import os
import re
import pandas as pd


def parse_filename(filename):
    base = os.path.basename(filename)
    if base.endswith(".results"):
        name = base[: -len(".results")]
    elif base.endswith(".fis"):
        name = base[: -len(".fis")]
    elif base.endswith(".txt"):
        name = base[: -len(".txt")]
        if name.startswith("allele_freqs_"):
            name = name[len("allele_freqs_") :]
    else:
        return None

    parts = name.split("_")
    if len(parts) % 2 != 0:
        return None

    params = {}
    for key, val in zip(parts[0::2], parts[1::2]):
        try:
            if val.isdigit():
                val_cast = int(val)
            else:
                val_cast = float(val)
        except ValueError:
            val_cast = val
        params[key] = val_cast

    return params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allele-freqs",
        required=True,
        nargs="+",
        help=("One or more paths to allele_freqs_*.txt files."),
    )
    parser.add_argument(
        "--output-pkl",
        required=True,
        help=("Path to output pickle file."),
    )
    args = parser.parse_args()

    rows = []

    for file in args.allele_freqs:
        params = parse_filename(file)
        if params is None:
            continue

        df_file = pd.read_csv(file)
        if not {"tick", "host_A", "symb_A"}.issubset(df_file.columns):
            continue

        ticks = df_file["tick"].values
        host_A_values = df_file["host_A"].values
        symb_A_values = df_file["symb_A"].values

        row = {
            "file": file,
            "ticks": ticks,
            "host_A_series": host_A_values,
            "symb_A_series": symb_A_values,
        }
        row.update(params)
        rows.append(row)

    if not rows:
        return

    allele_freqs_df = pd.DataFrame(rows)
    allele_freqs_df.to_pickle(args.output_pkl)


if __name__ == "__main__":
    main()