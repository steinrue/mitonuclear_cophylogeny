# Print the header to standard output (modify as needed to match your column names)
echo -e "chr\t$(head -1 ${1})"

# Loop through each file provided as input
for file in "$@"; do
    # Extract chromosome name based on the new file naming convention
    chromosome=$(basename "$file" | sed -E 's/.*chr([0-9XY]+).results/chr\1/')

    awk -v chr="$chromosome" '
        BEGIN { FS=OFS="\t" }
        NR > 1 { print chr, $0 }
    ' "$file"
done