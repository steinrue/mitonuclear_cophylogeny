import pandas as pd
import sys
import numpy as np

# Get the population sample list and mitochondrial file from input arguments
sample_list_path = sys.argv[1]
MITO_FILENAME = sys.argv[2]

with open(sample_list_path, 'r') as file:
    idsToExtract = [line.strip() for line in file]
data = pd.read_csv(MITO_FILENAME, sep="\t")

filtered_data = data[data["SampleID"].isin(idsToExtract)]

filtered_data = filtered_data.set_index("SampleID").reindex(idsToExtract).reset_index()

filtered_data["Range"] = filtered_data["Range"].str.replace(" ;", "")

# Remove the problematic variant "5899INS" from Found_Polys and Remaining_Polys
def clean_poly_string(poly_string):
    # Ensure the field is a string and remove the bad variant
    return " ".join([p for p in str(poly_string).split() if p != "5899INS!"])

filtered_data["Found_Polys"] = filtered_data["Found_Polys"].apply(clean_poly_string)
filtered_data["Remaining_Polys"] = filtered_data["Remaining_Polys"].apply(clean_poly_string)

filtered_data["Polys"] = filtered_data.apply(
    lambda row: '\t'.join(np.concatenate((row['Found_Polys'].split(), row['Remaining_Polys'].split()))),
    axis=1
)

required_cols = ["SampleID", "Range", "Haplogroup", "Polys"]
filtered_data = filtered_data[required_cols]

for row in filtered_data.values:
    print("\t".join(map(str, row)))
