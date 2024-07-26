```bash
# Find all sorted BAM files in the result_dorado directory
bam_file_sorted=$(find ./result_dorado/*_sorted.bam) # Path to sorted bam files

# Loop through each sorted BAM file found
for bam in $bam_file_sorted; do
    # Get the basename of the sorted BAM file, which is the patient ID
    id_patient=$(basename "$bam" _sorted.bam)

    # Execute the singularity command to run modkit pileup
    singularity exec --nv ./image_singularity/modkit.sif modkit pileup "$bam" ./result_modkit/"$id_patient.bed" --ref ./reference/*.fa --combine-strands --cpg
done
```
