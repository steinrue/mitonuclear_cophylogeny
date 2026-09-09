import sys
import pandas as pd
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <sample_list_path> <superpop_data>")
        sys.exit(1)

    sample_list_path = sys.argv[1]
    superpop_data = sys.argv[2]

    with open(sample_list_path, 'r') as file:
        idsToExtract = [line.strip() for line in file]

    df = pd.read_csv(superpop_data, sep=' ', na_values=[], keep_default_na=False)
    filtered_df = df[df['sample'].isin(idsToExtract)]
    filtered_df = filtered_df.set_index('sample').reindex(idsToExtract).reset_index()
    csv_output = filtered_df.to_csv(sep=' ', index=False)
    print(csv_output.strip())

if __name__ == "__main__":
    main()


