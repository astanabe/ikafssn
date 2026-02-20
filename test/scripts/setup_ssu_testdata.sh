#!/bin/bash
# Prepare derived test data from SSU_eukaryote_rRNA BLAST DB.
# Usage: ./setup_ssu_testdata.sh [SSU_PREFIX]
#
# Creates /tmp/ikafssn_ssu_test/ containing:
#   ssu_ambigdb.*     - BLAST DB with injected ambiguous bases
#   ssu_multivol_a.*  - BLAST DB subset (FJ876973.1 + 4 others)
#   ssu_multivol_b.*  - BLAST DB subset (GQ912721.1 + DQ235612.1 + 3 others)
#   queries.fasta     - 100bp query subsequences
#   seqidlist.txt     - Accession list for filter tests

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SSU_PREFIX="${1:-${SCRIPT_DIR}/../../db/SSU_eukaryote_rRNA}"
OUTDIR="/tmp/ikafssn_ssu_test"

# Check the source DB exists
if [ ! -f "${SSU_PREFIX}.nsq" ] && [ ! -f "${SSU_PREFIX}.00.nsq" ]; then
    echo "Error: SSU_eukaryote_rRNA not found at ${SSU_PREFIX}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

echo "=== Setting up SSU test data in ${OUTDIR} ==="

# ---- 1. Ambig DB ----
if [ ! -f "${OUTDIR}/ssu_ambigdb.nsq" ]; then
    echo "--- Building ssu_ambigdb ---"

    # Export all sequences to FASTA
    blastdbcmd -db "${SSU_PREFIX}" -dbtype nucl -entry all \
        -out "${OUTDIR}/SSU_eukaryote_rRNA.fasta"

    # Inject ambiguous bases via Python
    python3 -c "
import sys

infile = '${OUTDIR}/SSU_eukaryote_rRNA.fasta'
outfile = '${OUTDIR}/ssu_ambig.fasta'

# Read FASTA
records = {}
current_id = None
with open(infile) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            current_id = line.split()[0][1:]  # e.g. 'FJ876973.1'
            records[current_id] = {'header': line, 'seq': []}
        elif current_id:
            records[current_id]['seq'].append(line)

# Flatten sequences
for rid in records:
    records[rid]['seq'] = ''.join(records[rid]['seq'])

# Modify specific accessions
if 'FJ876973.1' in records:
    seq = list(records['FJ876973.1']['seq'])
    seq[100] = 'R'  # A|G at position 100
    records['FJ876973.1']['seq'] = ''.join(seq)

if 'GQ912721.1' in records:
    seq = list(records['GQ912721.1']['seq'])
    for i in range(50, 54):
        seq[i] = 'N'  # 4 Ns at positions 50-53
    records['GQ912721.1']['seq'] = ''.join(seq)

if 'DQ235612.1' in records:
    seq = list(records['DQ235612.1']['seq'])
    ambig = 'MRWSYKVHDB'
    for i, c in enumerate(ambig):
        seq[200 + i] = c  # 10 ambig bases at positions 200-209
    records['DQ235612.1']['seq'] = ''.join(seq)

# Write output
with open(outfile, 'w') as f:
    for rid in records:
        rec = records[rid]
        f.write(rec['header'] + '\n')
        seq = rec['seq']
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + '\n')
"

    makeblastdb -in "${OUTDIR}/ssu_ambig.fasta" -dbtype nucl \
        -out "${OUTDIR}/ssu_ambigdb" -parse_seqids

    # Clean up intermediate files
    rm -f "${OUTDIR}/SSU_eukaryote_rRNA.fasta" "${OUTDIR}/ssu_ambig.fasta"
    echo "  ssu_ambigdb created."
else
    echo "  ssu_ambigdb already exists, skipping."
fi

# ---- 2. Multi-volume DBs ----
if [ ! -f "${OUTDIR}/ssu_multivol_a.nsq" ]; then
    echo "--- Building ssu_multivol_a ---"
    blastdbcmd -db "${SSU_PREFIX}" -dbtype nucl \
        -entry "FJ876973.1,MH279652.1,MF325755.1,KY682881.1,KY682880.1" \
        -out "${OUTDIR}/ssu_multivol_a.fasta"
    makeblastdb -in "${OUTDIR}/ssu_multivol_a.fasta" -dbtype nucl \
        -out "${OUTDIR}/ssu_multivol_a" -parse_seqids
    rm -f "${OUTDIR}/ssu_multivol_a.fasta"
    echo "  ssu_multivol_a created."
else
    echo "  ssu_multivol_a already exists, skipping."
fi

if [ ! -f "${OUTDIR}/ssu_multivol_b.nsq" ]; then
    echo "--- Building ssu_multivol_b ---"
    blastdbcmd -db "${SSU_PREFIX}" -dbtype nucl \
        -entry "GQ912721.1,DQ235612.1,KY682879.1,KY682878.1,KY682877.1" \
        -out "${OUTDIR}/ssu_multivol_b.fasta"
    makeblastdb -in "${OUTDIR}/ssu_multivol_b.fasta" -dbtype nucl \
        -out "${OUTDIR}/ssu_multivol_b" -parse_seqids
    rm -f "${OUTDIR}/ssu_multivol_b.fasta"
    echo "  ssu_multivol_b created."
else
    echo "  ssu_multivol_b already exists, skipping."
fi

# ---- 3. Query FASTA ----
if [ ! -f "${OUTDIR}/queries.fasta" ]; then
    echo "--- Generating queries.fasta ---"
    {
        echo ">query_FJ876973 100bp from FJ876973.1 pos 101-200"
        blastdbcmd -db "${SSU_PREFIX}" -dbtype nucl \
            -entry "FJ876973.1" -range "101-200" -outfmt "%s"
        echo ">query_GQ912721 100bp from GQ912721.1 pos 51-150"
        blastdbcmd -db "${SSU_PREFIX}" -dbtype nucl \
            -entry "GQ912721.1" -range "51-150" -outfmt "%s"
    } > "${OUTDIR}/queries.fasta"
    echo "  queries.fasta created."
else
    echo "  queries.fasta already exists, skipping."
fi

# ---- 4. Seqidlist ----
if [ ! -f "${OUTDIR}/seqidlist.txt" ]; then
    echo "--- Generating seqidlist.txt ---"
    printf "FJ876973.1\nGQ912721.1\n" > "${OUTDIR}/seqidlist.txt"
    echo "  seqidlist.txt created."
else
    echo "  seqidlist.txt already exists, skipping."
fi

echo "=== Done. All test data in ${OUTDIR} ==="
