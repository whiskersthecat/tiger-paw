# -------------------------------
# 06_groupBlocks.py
# tiger-paw
# Peter Reifenstein
# 
# Use: group blocks by which primary haplotype classes they contain
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("blocks_file", help = "blocks file (.blocks)")

args = parser.parse_args()

ifile = open(args.blocks_file, 'r')
f_name = ".".join(args.blocks_file.split('.')[:-1])
ofile = open(f_name + ".grouped.blocks", 'w')
statsfile = open(f_name + ".stats.blocks", 'w')

grouped_reads = ({})

# -------------------------------
# 1. Read blocks

for full_line in ifile.readlines():
    line = full_line.split('\t')
    hap_seq = line[8].strip().split(' ')

    total_haps = 0

    groups = {}
    insdel_haps = {}

    for hap in hap_seq:
        group = hap[0]
        if group != "-" and group != "_":
            groups[group] = "present"
    id = ""
    for group in sorted(groups.keys()):
        id += group
    
    if len(id) == 0:
        id = "empty"

    if id not in grouped_reads:
         grouped_reads[id] = []
    
    grouped_reads[id].append(full_line)


sorted_ids = sorted(grouped_reads, key = lambda k: len(grouped_reads[k]), reverse = True)



# -------------------------------
# 2. Write blocks in order of which haplotype classes they contain

for id in sorted_ids:
    ofile.write("###\t" + id + "\n")
    ofile.write("##\tTotal " + str(len(grouped_reads[id])) + " reads \n")
    statsfile.write("###\t" + id + "\n")
    statsfile.write("##\tTotal " + str(len(grouped_reads[id])) + " reads \n\n")
    for line in grouped_reads[id]:
        ofile.write(line)
    ofile.write("\n")
