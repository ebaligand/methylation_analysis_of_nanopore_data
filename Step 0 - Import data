```bash
# Path to txt file
txt_file="/path/to/file/data_import.txt"

# Destination
destination="/path/to/destination/folder"

while IFS="     " read -r id_patient path; do
        mkdir $destination/$id_patient
        rsync $path/*.fast5 $destination/$id_patient --progress
done < $txt_file

# Replace "*.fast" xith "*.pod5" if your row files are in pod5 format 
```
