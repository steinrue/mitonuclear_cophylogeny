import pandas as pd
import subprocess
import sys
import os
from io import StringIO

pop_arg = sys.argv[1]
SAMPLE_FILENAME = sys.argv[2]
requested_pops = pop_arg.split("_")
idTable = pd.read_csv(SAMPLE_FILENAME, comment='#', sep="\t", header=None)

with open(SAMPLE_FILENAME, 'r') as file:
    for line in file:
        if line.startswith('#') and not line.startswith('##'):
            header = line.strip().lstrip('#').split('\t')
            break
idTable.columns = header

existing = set(idTable['POPULATION'].unique())
found_pops = [p for p in requested_pops if p in existing]
missing_pops = [p for p in requested_pops if p not in existing]

if missing_pops:
    sys.stderr.write(f"Warning: these populations were not found and will be skipped: {','.join(missing_pops)}\n")

if not found_pops:
    sys.exit("Error: none of the requested populations are present in the file.")

sampIdTable = idTable[idTable['POPULATION'].isin(found_pops)]

for sample in sampIdTable['SAMPLE_NAME'].str.strip():
    print(sample)

