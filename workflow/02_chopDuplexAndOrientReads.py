# -------------------------------
# 02_chopDuplexAndOrientReads.py
# tiger-paw
# Peter Reifenstein
# 
# Use: survey orientations of repeat segments and chop potential duplex reads
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("reads", help = "reads (.fa)")
parser.add_argument("blastn", help = "reformatted blastn results against reference (.sorted.blastN)")
parser.add_argument("--include_original", "-i", action = "store_true", help = "include original (unmodified) sequences in output fasta file")
args = parser.parse_args()

fasta_file = open(args.reads, 'r')
blastn_file = open(args.blastn, 'r')

output_fasta_file = open(args.reads + ".chop.oriented.fa", 'w')
output_readids_file = open(args.reads + ".readids", 'w')
output_report_file = open(args.reads + ".report", 'w')

include_original = args.include_original

RCmap = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
def reverse_complement(seq):
    rc = ""
    for n in seq[::-1]:
        rc += RCmap[n]
    return rc

raw_reads = {}

# -------------------------------
### 1. Read in fasta file into dictionary in memory

print("[process_reads] Reading reads into memory")
while True:
    name = fasta_file.readline().split()
    seq = fasta_file.readline().strip()
    if not seq:
        break
    name = name[0][1:]
    raw_reads[name] = seq

read_counts = [       0,         0,                0,                      0,                     0]
read_types  = ["forward", "reverse", "perfect_duplex", "bias_duplex_forward", "bias_duplex_reverse"]
read_abbrevs= ['f'      , 'r'      , 'dp'             ,  'dbf'                , 'dbr'                 ]

# -------------------------------
### 2. Look through the blastn file to analyze reads

print("[process_reads] Analyzing reads")

read_ctr = 0
previous_read_id = ""

orientationCTR = {"1/1" : 0, "1/-1" : 0}
orientationCoords = {"1/1" : 9999999999, "1/-1" : 9999999999}
currentOrientation = "0"

for result in blastn_file.readlines():

    result = result.split()
    read_id = result[1]

    if read_id != previous_read_id:
        read_ctr += 1
        if(read_ctr != 1 and (previous_read_id in raw_reads.keys())):
            seq = ""
            category = -1

            switch_point = max(orientationCoords["1/1"], orientationCoords["1/-1"])

            if(orientationCTR["1/1"] == 0 or orientationCTR["1/-1"] == 0):
                # 0, 1 Unidirectional read
                if(orientationCTR["1/-1"] == 0):
                    category = 0
                    seq = raw_reads[previous_read_id]

                if(orientationCTR["1/1"] == 0):
                    category = 1
                    seq = reverse_complement(raw_reads[previous_read_id])

            elif(orientationCTR["1/1"] == orientationCTR["1/-1"]):
                # 2 perfect duplex
                # ... take the first part
                if(currentOrientation == "1/1"):
                    seq = reverse_complement ( raw_reads[previous_read_id][:switch_point] )
                else:
                    seq = raw_reads[previous_read_id][:switch_point]

                category = 2

            elif(orientationCTR["1/1"] > orientationCTR["1/-1"]):
                # 3 bias duplex forward
                # ... cut out longer side
                if(currentOrientation == "1/1"):
                    seq = raw_reads[previous_read_id][switch_point:]
                else:
                    seq = raw_reads[previous_read_id][:switch_point]
                    
                category = 3

            elif(orientationCTR["1/1"] < orientationCTR["1/-1"]):
                # 4 bias duplex reverse
                # ... cut out longer side
                if(currentOrientation == "1/-1"):
                    seq = reverse_complement( raw_reads[previous_read_id][switch_point:] )
                else:
                    seq = reverse_complement( raw_reads[previous_read_id][:switch_point] )
                
                category = 4
            
            else:
                assert True, "error in logic"

            length = len(seq)

            output_fasta_file.write(">" + previous_read_id + "\tL:" + str(length) + "\n")
            output_fasta_file.write(seq + "\n")

            if(category >= 2 and include_original == True):
                length = len(raw_reads[previous_read_id])
                output_fasta_file.write(">" + previous_read_id + "_" + read_abbrevs[category] + "\tL:" + str(length) + "\n")
                output_fasta_file.write(raw_reads[previous_read_id] + "\n")

            output_readids_file.write(previous_read_id + '\t' + read_types[category] + '\n')
            read_counts[category] += 1

        orientationCTR = {"1/1" : 0, "1/-1" : 0}
        orientationCoords = {"1/1" : 9999999999, "1/-1" : 9999999999}
        currentOrientation = "0"
    
    previous_read_id = read_id

    start_coord = int(result[12])
    end_coord = int(result[13])
    orientation = result[9]

    try:
        orientationCTR[orientation] += 1
    
    except Exception as e:
        print("Exception: ", e)
        print("unknown orientation: " , orientation)
    
    if(currentOrientation == "0"):
        currentOrientation = orientation
    elif(currentOrientation == orientation):
        pass
    elif(currentOrientation != orientation):
        currentOrientation = orientation

    if (min(start_coord, end_coord) < orientationCoords[orientation]):
        orientationCoords[orientation] = min(start_coord, end_coord)

# -------------------------------
### 3. Output report

total = 0
for i in range(len(read_counts)):
    print(read_types[i], read_counts[i])
    output_report_file.write(read_types[i] + "\t" + str(read_counts[i]) + "\n")
    total += read_counts[i]

print("Total: " , total)
output_report_file.write("Total: " + "\t" + str(total) + "\n")

