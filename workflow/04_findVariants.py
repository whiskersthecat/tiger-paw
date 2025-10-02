# -------------------------------
# 04_findVariants.py
# tiger-paw
# Peter Reifenstein
# 
# Use: find variants in segments
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("sam_file", help = "segments aligned to reference (.sam)")
parser.add_argument("variants_file", help = "variants from highCoverageVariantCaller")
parser.add_argument("reference_sequence", help = "reference sequence (.fa)")

args = parser.parse_args()

sam_file = open(args.sam_file,'r')
variants_file = open(args.variants_file,'r')
ref_file = open(args.reference_sequence,'r')
group_name = args.sam_file

log_file = open(group_name + ".log.txt", "w")

print("[Find Variants] Group Name:", group_name)

line = ref_file.readline()
ref_seq = ref_file.readline().strip()



ref_len = len(ref_seq)
print("[INFO] reference length:", ref_len)

def resizeSTR(str, size):
    while(len(str) < size):
        str += "-"
    if(len(str) >= size):
        str = str[0:size]
    return str

read_names = ["Pos","Var","Len","Reference","Alternate","Frequency","Other_Count"]
alternate_count = {}
other_count = {}
Variations = {} #Variations[variant_type][location][variation]
VARIANT_TYPES = ["SNV","MNV","INS","DEL"]
for variant_type in VARIANT_TYPES:
    Variations[variant_type] = OrderedDict({})
    alternate_count[variant_type] = OrderedDict({})
    other_count[variant_type] = OrderedDict({})

# -------------------------------
## 1. Read in variants file

print("[Find Variants] Reading in Variants")
variants_file.readline() # consume header
for line in variants_file.readlines():
    line = line.split()
    var_type = line[1]
    location = int(line[0])
    variation = line[4]
    if(var_type == "DEL"):
        variation = line[3]
    
    if not (location in Variations[var_type]):
        Variations[var_type][location] = {}
        alternate_count[var_type][location] = {}
        other_count[var_type][location] = {}
    Variations[var_type][location][variation] = line[:5]
    alternate_count[var_type][location][variation] = 0
    other_count[var_type][location][variation] = 0

for variant_type in VARIANT_TYPES:
    print("[INFO]", len(Variations[variant_type]), "total", variant_type, "locations")
    log_file.write("[INFO]" + str(len(Variations[variant_type])) + "total" + variant_type + "locations\n")

# -------------------------------
## 2. Analyze reads in SAM file

print("[Find Variants] Parsing SAM CIGAR tokens for each read")

nreads = 0
reference_len = len(ref_seq)
for line in sam_file.readlines():

    line = line.split()
    if(line[0][0] == '@'): # skip header fields:
        continue

    name = line[0]
    flag = line[1]
    if (int(flag) & 2048 == 2048):
        print("[Find Variants] Skipping alignment with flag ", flag)
        log_file.write("[Find Variants] Skipping alignment with flag " + flag + "\n")
        continue
    cigar = line[5]
    seq = line[9]
    start_pos = int(line[3]) - 1

    insertions = {}
    deletions = {}

    # Generate the ALIGNED SEQUENCE
    align_seq = ""
    for i in range(0, int(start_pos)):   # if the alignment does not begin at position 0, buffer with - tokens:
        align_seq += "-"
    cur_num = ""
    pos = 0     # position in read sequence
    ref_pos = start_pos # position in reference sequence. should match with the length of align_seq
    
    if cigar == "*":
        print("[Find Variants] Skipping unaligned segment ", name)
        log_file.write("[Find Variants] Skipping unaligned segment " + name + "\n")
        continue

    for token in cigar:
        if token.isnumeric():
            cur_num += token
            continue
        try:
            num = int(cur_num)
        except:
            print("CIGAR PARSE ERROR:: UNKNOWN TOKEN:", token)
        new_pos = pos
        if token == "M" or token == "=" or token == "X":
            new_pos = pos + num
            align_seq += seq[pos:new_pos]
            ref_pos = ref_pos + num

        elif token == "I":
            new_pos = pos + num
            insertions[(ref_pos) + 1] = seq[pos:new_pos]

        elif token == "D":
            for i in range(0, num):
                align_seq += "-"
            deletions[(ref_pos) + 1] = ref_seq[ref_pos:num + ref_pos]
            ref_pos = ref_pos + num
        elif token == "S":
            new_pos = pos + num
        else:
            print("CIGAR PARSE ERROR:: UNKNOWN TOKEN:", token)
            exit()

        pos = new_pos
        cur_num = ""

    # if(nreads == 0):
    #     reference_len = len(align_seq)
    
    # BUFFER WITH - tokens if there is a soft clip at end of alignment
    while(len(align_seq) < reference_len):
        align_seq += "-"

    # QUALITY CHECK all reads should have the same align_seq length, which is the same as the reference
    if(len(align_seq) != reference_len):
        print("CRITICAL QUALITY CHECK FAILED:: aligned sequence not length of reference sequence for read", name, cigar, len(align_seq), reference_len)
        exit()

    read_names.append(name)

    for variant_type in VARIANT_TYPES:
        for location in Variations[variant_type]:
            for variation in Variations[variant_type][location]:

                variation_info = Variations[variant_type][location][variation]

                variation_len = len(variation)

                alternate = 4
                reference = 3
                if variant_type == "SNV":
                    compare_sequence = align_seq[location - 1]
                if variant_type == "MNV":
                    compare_sequence = align_seq[location - 1 : location - 1 + int(variation_info[2])]
                if variant_type == "INS":
                    if location in insertions:
                        compare_sequence = insertions[location]
                    else:
                        compare_sequence = variation_info[reference]
                if variant_type == "DEL":
                    alternate = 3
                    reference = 4
                    if location in deletions:
                        compare_sequence = deletions[location]
                    else:
                        compare_sequence = variation_info[reference]

                if compare_sequence == variation_info[alternate]:
                    # alternate found
                    Variations[variant_type][location][variation].append(variation_info[alternate])
                    alternate_count[variant_type][location][variation] += 1

                else:
                    # reference (or other)
                    if compare_sequence != variation_info[reference]:
                        other_count[variant_type][location][variation] += 1
                    
                    Variations[variant_type][location][variation].append(resizeSTR(compare_sequence, variation_len))
                   # Variations[variant_type][location][variation].append(variation_info[reference])

                    # Variations[variant_type][location][variation].append(variation_info[reference])
    
    nreads += 1

print("[INFO]", nreads, "total reads")
log_file.write("[INFO]" + str(nreads) + "total reads")


for variant_type in VARIANT_TYPES:
        for location in Variations[variant_type]:
            for variation in Variations[variant_type][location]:
                freq = 100 * alternate_count[variant_type][location][variation] / nreads
                Variations[variant_type][location][variation].insert(5, str(freq)[0:12])
                Variations[variant_type][location][variation].insert(6, other_count[variant_type][location][variation])

# -------------------------------
## 3. Generate large SNP table

print("[Find Variants] Generating Large Variant table")
largeOfile = open(group_name + ".variant_table.tab", "w")
for token in read_names:
    largeOfile.write(token + "\t")
largeOfile.write("\n")

for variant_type in VARIANT_TYPES:
    for location in Variations[variant_type].keys():
        for variation in Variations[variant_type][location]:
            for token in Variations[variant_type][location][variation]:
                largeOfile.write(str(token) + "\t")
            largeOfile.write("\n")

print("[Find Variants] Wrote all files")