# -------------------------------
# 07_stackedBlocksToSAM.py
# tiger-paw
# Peter Reifenstein
# 
# Use: fix the mapping coordinates of full length tandem repeat containing reads based on their manually aligned position during block stacking
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("stacked_blocks.stacks", help = "stacked blocks file, with proper headers")
parser.add_argument("tandem_alignment.sam", help = "reads aligned to tandem duplication of the reference in SAM format")
args = parser.parse_args()

labeled_file = open(args.manual_alignment.man, 'r')
sam_file = open(args.sam_alignment.sam, 'r')

MIN_SEGMENT_LENGTH = 0
REPEAT_LEN = 0
N_COMPONENTS = 0

contig = "NA"
contigs = {}
counts = {}

consensus_SAMs = {}
consensus_insertions = {}

def increment(d, k):
    if k not in d:
        d[k] = 1
    else:
        d[k] += 1

# -------------------------------
# 1. Parse the stacked blocks file and record the position and contig of each block

for line in labeled_file.readlines():
    if line == "\n":
        continue
    if line[0:2] == "# ":
        continue
    
    if line[0:4] == "####":
        contig = line[4:].strip()
        print("[CONTIG] " + contig)
        continue
    
    if line[0:8] == "X-DEFINE":
        arr = line.strip().split()
        consensus_name = arr[1]
        fpath = arr[2]
        try:
            consensus_sam_file = open(fpath, "r")
        except:
            print("[error] Could not open consensus sam file: ", fpath)
        seq = "@"
        while(seq[0] == "@"):
            seq = consensus_sam_file.readline().strip()
        consensus_SAMs[consensus_name] = seq.strip().split('\t')
        consensus_insertions[consensus_name] = {}
        print("[CONSENSUS] imported one consensus sequence named ", consensus_name)
        continue
    
    if line[0:12] == "X-REPEAT_LEN":
        arr = line.strip().split()
        REPEAT_LEN = int(arr[2])
        continue

    if line[0:12] == "X-COMPONENTS":
        arr = line.strip().split()
        N_COMPONENTS = len((arr[2].split(",")))
        continue

    if line[0:12] == "X-SEGMENTLEN":
        arr = line.strip().split()
        MIN_SEGMENT_LENGTH = int(arr[2])
        continue

    if (N_COMPONENTS == 0):
        print("[error] unknown number of components, add a special block at the top of the .stacks file:")
        print("X-COMPONENTS------------------------	-	telomere,copia_retrotransposon,right_linker")
        exit()

    if (MIN_SEGMENT_LENGTH == 0):
        print("[error] unknown minimum segment length, add a special block at the top of the .stacks file:")
        print("X-SEGMENTLEN------------------------	-	7000")
        exit()

    mode = "Align"
    if line[0:8] == "X-INSERT":
        mode = "Insert"

    arr = line.strip().split('\t')
    read_id = arr[0].strip()
    hap_seq = ""
    try:
        hap_seq = arr[5 + N_COMPONENTS]
    except:
        print("[error] Failed to properly extract haplotype list from block: ", arr)
    
    aln_shift = arr[4 + N_COMPONENTS]
    try:
        aln_shift = int(aln_shift)
    except:
        aln_shift = 0
    
    ## count number of spaces in front of the alignment
    spaces = 0
    while(hap_seq[spaces] == ' '):
        spaces += 1
    
    if (spaces % 3 != 0):
        print("[warning] block " + read_id + " has " + str(spaces) + " space characters, which is not a multiple of 3. Block position is rounded down.")
    
    position = int(spaces / 3)
    position += aln_shift

    if(mode == "Align"):
        if read_id in contigs:
            print("[warning] read " + read_id + " has already been assigned a position or contig")

        increment(counts, contig)
        contigs[read_id] = (contig, position)

    if(mode == "Insert"):
        type = arr[1]
        count = int(arr[2])
        depth = int(hap_seq[spaces:spaces+2])

        if not(type in consensus_SAMs):
            print("[error] consensus name ", type, " is not defined")
            exit()

        if not(contig in consensus_insertions[type]):
            consensus_insertions[type][contig] = []

        for n in range(count):
            # copy the consensus sequence alignment equal to the depth
            consensus_insertions[type][contig].append((position + n, depth))

if (REPEAT_LEN == 0):
    print("[error] repeat length is not defined, add a special block at the top of the .stacks file:")
    print("X-REPEAT_LEN------------------------	-	9847")
    exit()

CUTOFF_COORDINATE = REPEAT_LEN - MIN_SEGMENT_LENGTH

# -------------------------------
## 2. read the SAM file line by line and assign the read to an output file, adjusting the alignment start position

print("[output] making output files")

if not os.path.exists("drafts"):
    print("[output] making directory 'drafts'")
    os.mkdir("drafts")

o_sam_files = {}
for contig in counts:
    o_sam_files[contig] = open("drafts/" + contig + ".sam", 'w')
o_sam_files["unused"] = open("drafts/" + "unused" + ".sam", 'w')

report = open("drafts/report.txt", 'w')

for k in counts:
    print("[report] " + k + " has " + str(counts[k]) + " reads")
    report.write("[report] " + k + " has " + str(counts[k]) + " reads\n")

print("[output] reading SAM file")

unused = 0
counts = {}

adjusted = 0
unadjusted = 0

for line in sam_file.readlines():
    # write headers to every output file
    if(line[0] == "@"):
        for contig in o_sam_files:
            o_sam_files[contig].write(line)
        continue

    arr = line.strip().split('\t')
    read = arr[0].strip().replace("_","")
    position = int(arr[3])

    contig, shift = "", 0
    if read in contigs:
        contig, shift = contigs[read]
        del contigs[read]
    else:
        contig = "unused"
        unused += 1

    increment(counts, contig)
    
    if(contig != "unused"):
        left_justified_pos = (position % REPEAT_LEN)
        new_pos = left_justified_pos + (REPEAT_LEN * shift)
        # NEW --- check if there is a small segment in the front of the read that was removed:
        if (left_justified_pos > CUTOFF_COORDINATE):
            new_pos -= REPEAT_LEN
            adjusted += 1
        else:
            unadjusted += 1 

    arr[3] = str(new_pos)
    if(contig != "unused"):
        o_sam_files[contig].write('\t'.join(arr) + '\n')

print("[report] total " + str(unused) + " unused reads")
report.write("[report] total " + str(unused) + " unused reads\n")

# print("[report] total " + str(adjusted) + " adjusted reads")
# report.write("[report] total " + str(adjusted) + " adjusted reads\n")

# print("[report] total " + str(unadjusted) + " unadjusted reads")
# report.write("[report] total " + str(unadjusted) + " unadjusted reads\n")


for k in counts:
    print("[report] " + k + " has " + str(counts[k]) + " reads")
    report.write("[report] " + k + " has " + str(counts[k]) + " reads\n")


# -------------------------------
# 3. write the consensus insertions

print("[output] writing consensus insertions")

for consensus in consensus_insertions:
    seg_count = 0
    for contig in consensus_insertions[consensus]:
        for insertion in consensus_insertions[consensus][contig]:
            position, depth = insertion
            seg_count += 1
            depth_count = 0
            for i in range (depth):

                ### finish this part
                arr = consensus_SAMs[consensus].copy()
                depth_count += 1

                read = arr[0]

                new_pos = int(arr[3]) + (REPEAT_LEN * position)
                arr[3] = str(new_pos)
                new_id = read + "_" + str(seg_count) + "_" + str(depth_count)
                arr[0] = str(new_id)

                o_sam_files[contig].write('\t'.join(arr) + '\n')

    print("wrote a total of ", seg_count, "alignments for consensus ", consensus)


print("[output] done writing SAM file")

if(len(contigs.keys()) > 0):
    print("[warning] these reads were in the stacked blocks but were not found in SAM file:")
    print(contigs)
