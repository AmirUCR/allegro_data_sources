#!/bin/sh
DATABASES_DIR="inputs/databases/"
REF_SPECIES_PATH="inputs/reference/"
REF_SPECIES_PATH="${REF_SPECIES_PATH}$(ls -t $REF_SPECIES_PATH | head -n 1)"
OUTPUT_DIR="outputs/"

echo "Working..."
for DB in $(ls $DATABASES_DIR)
do
    SPECIES_NAME="$(echo "$DB" | cut -d'.' -f1)"
    DB_NAME="${DATABASES_DIR}${SPECIES_NAME}"
    OUT="${OUTPUT_DIR}${SPECIES_NAME}.tsv"

    diamond blastp -q $REF_SPECIES_PATH -d $DB_NAME -o $OUT --very-sensitive --max-target-seqs 1 --quiet
done
echo "Done."