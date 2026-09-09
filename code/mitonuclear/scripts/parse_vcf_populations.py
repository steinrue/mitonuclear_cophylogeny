import pandas as pd
import subprocess
import sys
import os

BCFTOOLS_BIN = 'bcftools'

sample_list_path = sys.argv[1]
vcfFile = sys.argv[2]
outFile = sys.argv[3]
with open(sample_list_path, 'r') as file:
    idsToExtract = list([line.strip() for line in file])

bcfToolsCmd = f'{BCFTOOLS_BIN} view --force-samples -s {",".join(idsToExtract)} {vcfFile} | {BCFTOOLS_BIN} annotate -x INFO - > {outFile}'

print(bcfToolsCmd)
subprocess.run(bcfToolsCmd, shell=True)