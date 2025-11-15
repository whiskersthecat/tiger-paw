# -------------------------------
# 03_extractSegments.py
# tiger-paw
# Peter Reifenstein
# 
# Use: extract repeat segments from input fasta file and blastN results, including large mystery sequences (undetected segments) between blast hits
# Note: in the genome assembly of the NOR region in Lettuce, these undetected segments are found to be Copia retrotransposons
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("reads", help = "reads (.fa)")
parser.add_argument("blastn", help = "reformatted blastn results against reference (.sorted.blastN)")
parser.add_argument("--edgeTerminalMinimumLength", "-l", help = "minimum length to include terminal segment in output (int)", default = 1000)
parser.add_argument("--undetectedSegmentMinimumLength", "-m", help = "minimum length to include undetected segment between two segments (int)", default = 100)
parser.add_argument("--strictFirstSegment", action = "store_true", help = "remove first segment if alignment coordinate starts at 1 (even if its large)")

args = parser.parse_args()

fasta_file = open(args.reads, 'r')
blastn_file = open(args.blastn, 'r')

segment_file = open(args.reads + ".segmented.fa", 'w')
undetected_file = open(args.reads + ".undetected.fa", 'w')
all_file = open(args.reads + ".allregions.fa", 'w')
log_file = open(args.reads + ".report", 'w')


MIN_TERMINAL_EDGE_LEN = int(args.edgeTerminalMinimumLength)
MIN_UNDETECTED = int(args.undetectedSegmentMinimumLength)

STRICT_FIRST_SEGMENT = bool(args.strictFirstSegment)

print(STRICT_FIRST_SEGMENT)

RCmap = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def reverse_complement(seq):
    rc = ""
    for n in seq[::-1]:
        rc += RCmap[n]
    return rc

def flush_segments(current_segments, read_id):
    skipped_edges = 0
    num_segs = len(current_segments)
    global total_seg, total_undetected
    for n, segmentinfo in enumerate(current_segments):
        segment_id, type, segment_len, start_coord, end_coord = segmentinfo

        if n == 0 or n == num_segs - 1:
            if n == num_segs - 1 and segment_len < MIN_TERMINAL_EDGE_LEN:
            
                skipped_edges += 1
                continue

            if n == 0 and STRICT_FIRST_SEGMENT and start_coord == 1:
                # print("Skipped with start coord == 1")
                skipped_edges += 1
                continue
            ### used for hifi reads
            if n == 0 and (segment_len < MIN_TERMINAL_EDGE_LEN):
                skipped_edges += 1
                continue

            ### used for other reads
            # if n == 0 and segment_len < MIN_TERMINAL_EDGE_LEN:
                
        try:
            segment_seq = raw_reads[read_id][start_coord : end_coord]
        except:
            print("Warning: read_id ", read_id, "not found in fa file")
            return skipped_edges

        if(segment_len < 0):
            segment_seq =  reverse_complement(raw_reads[read_id][end_coord : start_coord])
            # print("reverse complement : ", segment_seq)

        if(type == "U"):
            undetected_file.write(">" + segment_id + "\n" + segment_seq + "\n")
            all_file.write(">" + segment_id + "\tC:UNDETECTED" + "\n" + segment_seq + "\n")
            total_undetected += 1

        if(type == "R"):
            segment_file.write(">" + segment_id + "\n" + segment_seq + "\n")
            all_file.write(">" + segment_id + "\tC:RIBO" + "\n" + segment_seq + "\n")
            total_seg += 1

    return skipped_edges

raw_reads = {}


# -------------------------------
# 1. read in fasta file into dictionary in memory

print("[extract_segments] Reading reads into memory")
while True:
    name = fasta_file.readline().split()
    seq = fasta_file.readline().strip()
    if not seq:
        break
    name = name[0][1:]
    raw_reads[name] = seq


# -------------------------------
# 2. extract segments from fasta file

print("[extract_segments] Extracting segments")
read_ctr = 0
previous_read_id = ""
total_undetected = 0
total_seg = 0

skipped_edges = 0

current_segments = []

for result in blastn_file.readlines():
    result = result.split()
    read_id = result[1]
    if read_id != previous_read_id:
        skipped_edges += flush_segments (current_segments, previous_read_id)
        read_ctr += 1
        segment_ctr = 0
        previous_end_coord = 0
        current_segments = []
    
    previous_read_id = read_id
    segment_ctr += 1

    start_coord = int(result[12])
    end_coord = int(result[13])

    segment_len = min(start_coord, end_coord) - previous_end_coord
    if(start_coord > end_coord):
        segment_len = end_coord - previous_end_coord

    if segment_ctr > 1 and abs(segment_len) >= MIN_UNDETECTED: # undetected segment between NOR segments

        # note: the coordinates of previous undetected segments are wrong

        segment_id = read_id + "_" + str(read_ctr) + "_" + str(segment_ctr) + "\tS:" + str(start_coord) + "\tE:" + str(previous_end_coord) + "\tL:" + str(segment_len)
        current_segments.append((segment_id, "U", segment_len, previous_end_coord, start_coord))

        # segment_seq = raw_reads[read_id][previous_end_coord : start_coord]
        # undetected_file.write(">" + segment_id + "\n" + segment_seq + "\n")
        # all_file.write(">" + segment_id + "\tC:UNDETECTED" + "\n" + segment_seq + "\n")
        segment_ctr += 1
        # total_undetected += 1

    # total_seg += 1
    segment_len = end_coord - start_coord
    
    previous_end_coord = max(start_coord, end_coord)
    segment_id = read_id + "_" + str(read_ctr) + "_" + str(segment_ctr) + "\tS:" + str(start_coord) + "\tE:" + str(end_coord) + "\tL:" + str(segment_len)
    
    current_segments.append((segment_id, "R", segment_len, start_coord, end_coord))

    # segment_seq = raw_reads[read_id][start_coord : end_coord]
    # segment_file.write(">" + segment_id + "\n" + segment_seq + "\n")
    # all_file.write(">" + segment_id + "\tC:RIBO" + "\n" + segment_seq + "\n")

# -------------------------------
# 3. output summary

log_file.write("Total read       :   " + str(read_ctr) + "\n")
log_file.write("Total segments   :   " + str(total_seg) + "\n")
log_file.write("Total undetected :   " + str(total_undetected) + "\n")
log_file.write("Total skipped edges: " + str(skipped_edges) + "\n")

print("Total read       :", read_ctr)
print("Total segments   :", total_seg)
print("Total undetected :", total_undetected)
print("Total skipped edges:", skipped_edges)

