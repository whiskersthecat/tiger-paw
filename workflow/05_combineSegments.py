# -------------------------------
# 05_combineSegments.py
# tiger-paw
# Peter Reifenstein
# 
# Use: combine classified segments into the reads they came from
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("classified_file", help = "list of segment ids and their haplotypes (.Traces.Classified.tsv)")

args = parser.parse_args()

class_file = open(args.classified_file, 'r')
reconstructed_file = open(args.classified_file + ".haplotype.list.tab", 'w')

MISSING_SEG = "__"

cur_read_name = ""
prev_read_name = ""
hap_seq = ""
read_index = 0

print("[INFO] reading segment classifications")

for line in class_file.readlines():
    line = line.split()
    read = line[0].split('_')
    cur_read_name = '_'.join(read[:-2])
    if cur_read_name != prev_read_name and prev_read_name != "":
        reconstructed_file.write(prev_read_name + "\t" + hap_seq + "\n")
        hap_seq = ""
        read_index = 0
    read_index += 1
    while read_index < int(read[-1]):
        read_index += 1
        ### NEW: doesn't add missing haplotype for first segment
        if(read_index != 2):
            hap_seq += MISSING_SEG + " "
    
    hap_seq += line[-2] + " "
    prev_read_name = cur_read_name

reconstructed_file.write(prev_read_name + "\t" + hap_seq + "\n")

print("[INFO] done")
