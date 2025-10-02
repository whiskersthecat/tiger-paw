
MINIMUM_N_SEGMENTS=3
MOST_SEGMENTS=25

spacing=$((MOST_SEGMENTS*3))
min_spaces=$(((MINIMUM_N_SEGMENTS-1)*3))


seperator_index=-1
args=("$@")
for i in "${!args[@]}"; do
  if [[ "${args[i]}" == "--" ]]; then
    sep_index=$i
  fi
done

if (( sep_index < 0 )); then
  echo "Usage: $0 dataset_name <list of blast coverage results .coverage> -- <list of combined haplotype files .combined>"
  echo "E.g.: $0 short_1 <reference.vs....bed.coverage component1.vs....bed.coverage component2.vs....bed.coverage ...> -- <var1....combined var2....combined "
  exit 1
fi

dataset_name=$1
dataset_abbrev=${dataset_name:0:1}
echo "[input] dataset name: $dataset_name"
echo "[input] dataset abbreviation name: ${dataset_abbrev}"

reference_coverage_file=$2
echo "[input] reference coverage file: $reference_coverage_file"

coverage_files=("${args[@]:2:sep_index-2}")
# echo "[input] Component coverage files: ${coverage_files[*]}"

echo "[output] making directory "blocks""
mkdir -p blocks


echo "[part 1] cutting from reference hits coverage table"
cut -f 1,3,7 ${reference_coverage_file} | sort > blocks/${dataset_name}.blocks


for coverage_file in "${coverage_files[@]}"; do
  echo "[part 1] processing coverage file $coverage_file"
  cut -f 1,7 ${coverage_file} | sort > ${coverage_file}.tmp
  join -t $'\t' blocks/${dataset_name}.blocks ${coverage_file}.tmp > blocks/${dataset_name}.blocks.tmp
  mv blocks/${dataset_name}.blocks.tmp blocks/${dataset_name}.blocks
done

echo "[part 1] reformatting intermediate file with all coverages"

# cp blocks/${dataset_name}.blocks blocks/${dataset_name}.blocks

awk -v VAR="$dataset_abbrev" 'BEGIN {FS=OFS="\t"} {
  printf ("%s\t%s\t%-7s\t", $1, VAR, $2); other = 1; fields = NF;
  for(i = 3; i <= fields; i++) { other -= $i; if(i > 3) {if ($i > 0.001) {printf("\t%02.1i", $i * 100)} else {printf("\t--")};} };
  if (other > 0.01) {printf("\t%02.1i", other * 100)} else {printf("\t--")}; 
  printf("\n")} ' blocks/${dataset_name}.blocks > blocks/${dataset_name}.blocks.tmp
mv blocks/${dataset_name}.blocks.tmp blocks/${dataset_name}.blocks


combined_files=("${args[@]:sep_index+1}")
# echo "[part 2] Combined files: ${combined_files[*]}"


echo "[part 2] setting max number of haplotypes to $MOST_SEGMENTS"


echo "[part 2] setting spacing $spacing"

for combined_file in "${combined_files[@]}"; do
  echo "[part 2] processing combined file $combined_file"
  sort ${combined_file} | awk -v VAR="$min_spacer" 'BEGIN {FS=OFS="\t"} length($2) > VAR' | awk -v VAR="$spacing" 'BEGIN {FS=OFS="\t"} {printf("%s\t%-*s\n",$1,VAR,$2)}'   > ${combined_file}.sorted 
  join -t $'\t' blocks/${dataset_name}.blocks ${combined_file}.sorted > blocks/${dataset_name}.blocks.tmp
  mv blocks/${dataset_name}.blocks.tmp blocks/${dataset_name}.blocks
done

echo "[formatting] doing final reformatting"

sort blocks/${dataset_name}.blocks -k3,3nr | awk 'BEGIN {FS=OFS="\t"} {$1 = sprintf("%-40s", $1); print}'> blocks/${dataset_name}.blocks.tmp

mv blocks/${dataset_name}.blocks.tmp blocks/${dataset_name}.blocks

