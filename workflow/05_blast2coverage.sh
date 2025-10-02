
blast_file="$1"
targets_fa="$2"

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <special_arg> [file1 file2 ...]" >&2
  exit 1
fi

if ! command -v bedtools &>/dev/null; then
  echo "Error: bedtools not found" >&2
  exit 1
fi

# 1. write the bed file of the reads

if [ ! -f "$targets_fa.bed" ]; then
  echo "Creating bed files from reads"
  awk 'BEGIN {FS=OFS="\t"} {printf("%s\t", substr($1, 2)); getline; printf("0\t%d\n", length($1))}' "$targets_fa" > "$targets_fa.bed"
fi

# 2. write the bed file of the hits
echo "Creating bed files from blastn file"

awk 'BEGIN {FS=OFS="\t"} !/^#/ {s = $2; a = $13; b = $14; if (a <= b) {start0 = a - 1; end0   = b} else {start0 = b - 1; end0 = a;} print s"\t"start0"\t"end0 }' "$blast_file" > "${blast_file}.hits.bed"

# 3. sort bed files
echo "Sorting bed files"

sort -k1,1 -k2,2n "${targets_fa}.bed" > "${targets_fa}.sorted.bed"
sort -k1,1 -k2,2n "${blast_file}.hits.bed" > "${blast_file}.hits.sorted.bed"

# 4. merge them with bedtools
echo "Calculating coverage in bed files with bedtools"

bedtools coverage -a "${targets_fa}.sorted.bed" -b "${blast_file}.hits.sorted.bed" | sort -k7,7nr > "${blast_file}.hits.bed.coverage"

