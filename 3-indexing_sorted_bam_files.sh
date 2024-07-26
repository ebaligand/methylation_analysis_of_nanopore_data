```bash
# Find all BAM files in the result_dorado directory
bam_files=$(find ./result_dorado -name "*.bam") # Path to BAM files after basecalling

# Loop through each BAM file found
for bam in $bam_files; do
    # Get the basename of the BAM file, which is the patient ID
    id_patient=$(basename "$bam" .bam)
    # Define the sorted BAM file name
    sorted_bam="./result_dorado/${id_patient}_sorted.bam"

    # Execute the singularity command to sort the BAM file
    singularity exec --nv ./image_singularity/samtools.sif samtools sort -o "$sorted_bam" "$bam"

    # Execute the singularity command to index the sorted BAM file
    singularity exec --nv ./image_singularity/samtools.sif samtools index "$sorted_bam"
done

echo "Sorting and indexing of all BAM files completed."
```
