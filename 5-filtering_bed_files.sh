```bash
# Find all BED files in the result_modkit directory
bed_files=$(find ./result_modkit/*.bed) # Path to bed files

# Loop through each BED file found
for bed in $bed_files; do
    # Get the basename of the BED file, which is the patient ID
    id_patient=$(basename "$bed" .bed)

    # Use awk to filter rows where the fourth column is "m" and save to a new filtered BED file
    awk '$4=="m"' "$bed" > "./result_modkit/${id_patient}_filtered.bed"

    # Get the new filtered BED file path
    filtered_bed="./result_modkit/${id_patient}_filtered.bed"

    # Use awk to print the patient ID, chromosome and position, and the 11th column to a new text file
    awk -v pid="$id_patient" '
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|[XYM])$/ {
            print pid, $1":"$2, $11
        }
    ' "$filtered_bed" > "./result_modkit/${id_patient}_modified.txt"
done
# You can eliminate chromosomes from the analysis by deleting them on line 18 (exemple to exclude the Y chromosome : chr([1-9]|1[0-9]|2[0-2]|[XM])
# To keep the hydroxymethylation information, replace this: '$4=="m"' with that: '$4=="m" || $4=="h"'
```
