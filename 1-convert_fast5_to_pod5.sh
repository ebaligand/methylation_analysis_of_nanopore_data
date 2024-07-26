```bash
# Find all directories at the specified path
data_folder=$(find /path/to/data_fast5 -mindepth 1 -maxdepth 1 -type d) # Replace /path/to/data_fast5 with your path foward fast5 files
destination="./data_pod5"

# Loop through each folder found
for folder in $data_folder; do
    # Get the basename of the folder, which is the patient ID
    id_patient=$(basename "$folder")
    # Create a directory for each patient ID in the destination directory
    mkdir -p "$destination/$id_patient"

    # Execute the singularity command to convert fast5 files to pod5 format
    singularity exec --nv ./image_singularity/fast5_to_pod5.sif pod5 convert fast5 "$folder"/*.fast5 --output "$destination/$id_patient" --one-to-one "$folder"
done
```
