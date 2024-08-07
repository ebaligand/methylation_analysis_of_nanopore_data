```bash
# Find all directories in the data_pod5 directory
data_folder=$(find ./data_pod5 -mindepth 1 -maxdepth 1 -type d) # ./data_pod5 is the path to pod5 files

# Loop through each folder found
for folder in $data_folder; do
    # Get the basename of the folder, which is the patient ID
    id_patient=$(basename "$folder")

    # Execute the singularity command to run the dorado basecaller
    singularity exec --nv ./image_singularity/dorado.sif dorado basecaller --reference ./reference/*.fa /models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 "$folder" --modified-bases-models /models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v1 > ./result_dorado/"$id_patient.bam"
done
```
