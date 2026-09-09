import sys
if len(sys.argv) != 2:
    print("Usage: python script.py input_file")
    sys.exit(1)

input_file = sys.argv[1]

with open(input_file, 'r') as file:
    lines = file.readlines()
lines = lines[1:]

new_lines = []
for line in lines:
    columns = line.strip().split()
    columns[4] = "NA"
    rearranged_columns = [columns[1], columns[5], columns[6], columns[4]]
    new_line = ' '.join(rearranged_columns)
    new_lines.append(new_line)
header = "sample population group sex"
print(header)
print('\n'.join(new_lines))
