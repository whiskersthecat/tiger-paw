# -------------------------------
# 10_computeKmers.sh
# tiger-paw
# Peter Reifenstein
# 
# Use: determine fraction of kmers in reads in each genome assembly
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

# INPUT
genomes=( genome/genome.fa genome/noTRC.fa genome/Contig1.fa genome/originalContig1.fa )
datasets=( reads/long.fa reads/short.fa )

# PARAMETERS
k=25
threads=40
mincountreads=2

echo "[compute_kmers] threads=${threads}"
echo "[compute_kmers] k=${k}"
echo "[compute_kmers] mincountreads=${mincountreads}"
echo ""

weighted_output="kmers_weighted_m${mincountreads}_k${k}.tsv"
unique_output="kmers_unique_m${mincountreads}_k${k}.tsv"

echo "[compute_kmers] output file: ${weighted_output}"
echo "[compute_kmers] output file: ${unique_output}"
echo ""

printf "_WEIGHTED_\t" > $weighted_output
printf "__UNIQUE__\t" > $unique_output
for genome in "${genomes[@]}"; do
  printf "${genome}\t" >> $weighted_output
  printf "${genome}\t" >> $unique_output
done 
printf "\n" >> $weighted_output
printf "\n" >> $unique_output

for dataset in "${datasets[@]}"; do
  printf "${dataset}\t" >> $weighted_output
  printf "${dataset}\t" >> $unique_output
  for genome in "${genomes[@]}"; do
    echo "[dataset] $dataset [genome] $genome"
    genome_name="${genome}_k${k}"
    dataset_name="${dataset}_m${mincountreads}_k${k}"
    echo "[0] output genome name: $genome_name"
    echo "[0] output dataset name: $dataset_name"

    if [ ! -f "${dataset_name}.kmc_pre" ]; then
      echo "[0] counting $dataset"
      kmc -$threads -ci${mincountreads} -k$k -cs16777216 -fm "$dataset.fa" "${dataset_name}" "tmp/" > "${dataset_name}.log"
    fi
    if [ ! -f "$genome_name.kmc_pre" ]; then
      echo "[0] counting $genome"
      kmc -t10 -ci1 -b -k$k -cs65536 -fm "$genome.fa" "$genome_name" "tmp/" > "${genome_name}.log"
    fi

    if [ ! -f "$genome.stats" ]; then
      echo "[0] computing assembly stats for $genome"
      assembly-stats "$genome.fa" > "$genome.stats"
    fi
    
    if [ ! -f "$dataset.stats" ]; then
      echo "[0] computing assembly stats for $dataset"
      assembly-stats "$dataset.fa" > "$dataset.stats"
    fi

    if [ ! -f "${dataset_name}.intersect.$genome_name.kmc_pre" ]; then
      echo "[1] intersecting"
      kmc_tools simple "${dataset_name}" "$genome_name" intersect "${dataset_name}.intersect.$genome_name"
      echo ""
      kmc_dump "${dataset_name}.intersect.$genome_name" "${dataset_name}.intersect.$genome_name.txt"
    fi
    if [ ! -f "${dataset_name}.txt" ]; then
      echo "[2] dumping"
      kmc_dump "${dataset_name}" "${dataset_name}.txt"
    fi
    # if [ ! -f "$genome_name.txt" ]; then
    #   echo "[2] dumping"
    #   kmc_dump "$genome_name" "$genome_name.txt"
    # fi

    echo "[3] counting"

    if [ ! -f "${dataset_name}.sum" ]; then
      awk 'BEGIN {sum = 0} {sum += $2} END {print(sum)}' "${dataset_name}.txt" > "${dataset_name}.sum"
    fi
    if [ ! -f "${dataset_name}.intersect.${genome_name}.sum" ]; then
      awk 'BEGIN {sum = 0} {sum += $2} END {print(sum)}' "${dataset_name}.intersect.${genome_name}.txt" > "${dataset_name}.intersect.$genome_name.sum"
    fi

    if [ ! -f "${dataset_name}.nlines" ]; then
      wc -l "${dataset_name}.txt" | awk '{print $1}' > "${dataset_name}.nlines"
    fi
    if [ ! -f "${dataset_name}.intersect.${genome_name}.nlines" ]; then
      wc -l "${dataset_name}.intersect.${genome_name}.txt" | awk '{print $1}' > "${dataset_name}.intersect.${genome_name}.nlines"
    fi

    intersect_count=$(cat "${dataset_name}.intersect.$genome_name.sum")
    dataset_count=$(cat "${dataset_name}.sum")
    # echo "intersect_count: ${intersect_count}"
    # echo "dataset_count: ${dataset_count}"
    if [ "$dataset_count" -gt 0 ]; then
      percent=$( echo "100 * $intersect_count / $dataset_count" | bc -l )
    else
      percent="NA"
    fi
    echo "$percent" > "${dataset_name}.intersect.$genome_name.weighted.percent"
    printf "%.10f\t" ${percent} >> $weighted_output

    intersect_count=$(cat "${dataset_name}.intersect.$genome_name.nlines")
    dataset_count=$(cat "${dataset_name}.nlines")
    if [ "$dataset_count" -gt 0 ]; then
      percent=$( echo "100 * $intersect_count / $dataset_count" | bc -l )
    else
      percent="NA"
    fi
    echo "$percent" > "${dataset_name}.intersect.$genome_name.unique.percent"
    printf "%.10f\t" ${percent} >> $unique_output
    echo ""
  done
  printf "\n" >> $weighted_output
  printf "\n" >> $unique_output
done




