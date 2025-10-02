# -------------------------------
# 05_classifySegmentHaplotypes.py
# tiger-paw
# Peter Reifenstein
# 
# Use: classify segments by their closest matching haplotype
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

from collections import OrderedDict
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("traces_file", help = "list of segment ids and their traces (.Traces)")
parser.add_argument("haplotypes_file", help = "list of haplotypes (ID, sequence) (.tsv)")

args = parser.parse_args()


read_variants = open(args.traces_file,'r')
reference_variants = open(args.haplotypes_file,'r')

output_file = open(args.traces_file + ".Classified.tsv", 'w')

stats_file = open(args.traces_file + ".STATS.tsv", 'w')
stats_file.write("Haplotype\tCount\tPercentage\n")

# 1. Load dictionary queried by haplotype into memory
# -------------------------------

print("[INFO] Loading haplotype traces into memory")

reference_haplotype_names = OrderedDict({})     # trace      -> haplotype
reference_haplotype_count = OrderedDict({})     # haplotype  -> count
haplotype_name_len = 0
for line in reference_variants.readlines():
    line = line.split()
    reference_haplotype_names[line[1]] = line[0]
    reference_haplotype_count[line[0]] = 0
    reference_haplotype_count[line[0].lower()] = 0
    haplotype_name_len = len(line[0])

no_haplotype_name = "-"
while len(no_haplotype_name) < haplotype_name_len:
    no_haplotype_name += "-"
reference_haplotype_count[no_haplotype_name] = 0


# 2. Categorize reads
# -------------------------------

print("[INFO] Categorizing segments")

total_reads = 0
for line in read_variants.readlines():
    if(total_reads % 1000 == 0):
        print("categorized ", total_reads, " reads")
    
    line = line.split()
    # default is uncategorized
    read_haplotype = no_haplotype_name
    read_trace = line[1]
    read_trace = read_trace.replace("N", ".")

    read_id = line[0]

    # 2A. PERFECT MATCH
    # check for perfect haplotype categorization
    for ref_trace in reference_haplotype_names:
        if re.fullmatch(read_trace, ref_trace):
            read_haplotype = reference_haplotype_names[ref_trace]        

    # 2B. CLOSEST MATCH
    # count how many mismatches each haplotype has to this one. Assign to haplotype with fewest mismatches.
    fewest_mismatches = 9999
    best_haplotype = "NA"
    if read_haplotype == no_haplotype_name:
        for ref_trace in reference_haplotype_names:
            mismatch = 0
            for (i, nucleotide) in enumerate(ref_trace):
                if (read_trace[i] != nucleotide):
                    mismatch += 1
            
            if (mismatch < fewest_mismatches):
                best_haplotype = reference_haplotype_names[ref_trace].lower()
                fewest_mismatches = mismatch
         
        read_haplotype = best_haplotype
    
    output_file.write(line[0] + '\t' + read_haplotype  + '\t' + read_trace.strip() + '\n')
    reference_haplotype_count[read_haplotype] += 1
    total_reads += 1

# 3. Write stats
# -------------------------------

print("[INFO] Writing stats")

all_haplotype_names = reference_haplotype_names.copy()

for haplotype in reference_haplotype_count:
    count = reference_haplotype_count[haplotype]
    percentage = 100 * count / total_reads
    stats_file.write(haplotype + "\t" + str(count) + "\t" + str(percentage)[0:10] + '\n')

stats_file.write("TOTAL\t" + str(total_reads) + "\t100\n")

print("[INFO] Done")