mkdir -p gtdb_r226 && cd gtdb_r226

# Proteomes (amino acid sequences) from the MAGs
curl -O https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/genomic_files_reps/gtdb_proteins_aa_reps_r226.tar.gz

# Metadata for bac120 genomes
curl -O https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/bac120_metadata_r226.tsv.gz

# Check file sizes and compression format
ls -lh
file *.gz

awk -F'\t' '
BEGIN{ IGNORECASE=1; OFS="\t" }
NR==1{
    # Header row
    print $0 > "soil_metadata.tsv"
    next
}
# Filter: column 62 (ncbi_isolation_source) contains "soil"
# and column 63 (ncbi_lat_lon) is valid
$62 ~ /soil/ && $63 != "" && tolower($63) != "none" {
    # Normalize accession ID
    id=$1
    gsub(/^(GB_|RS_)/, "", id)
    sub(/\.[0-9]+$/, "", id)
    $1=id

    # Include all metadata
    print >> "soil_metadata.tsv"

    # Collect accessions
    print $1 >> "soil_ids.txt"
}
' bac120_metadata_r226.tsv

# Count selected entries (7,341)
tail -n +2 soil_metadata.tsv | wc -l

