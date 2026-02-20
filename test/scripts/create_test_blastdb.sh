#!/bin/bash
# Create test BLAST DB from small.fasta
# Usage: ./create_test_blastdb.sh [output_dir]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTDATA_DIR="${SCRIPT_DIR}/../testdata"
OUTPUT_DIR="${1:-${TESTDATA_DIR}}"

FASTA="${TESTDATA_DIR}/small.fasta"
DB_PREFIX="${OUTPUT_DIR}/testdb"

if [ ! -f "${FASTA}" ]; then
    echo "Error: ${FASTA} not found"
    exit 1
fi

echo "Creating test BLAST DB from ${FASTA}..."
makeblastdb -in "${FASTA}" -dbtype nucl -out "${DB_PREFIX}" -parse_seqids

if [ $? -eq 0 ]; then
    echo "Test BLAST DB created at ${DB_PREFIX}"
    echo "Files:"
    ls -la "${DB_PREFIX}".*
else
    echo "Error: makeblastdb failed"
    exit 1
fi
