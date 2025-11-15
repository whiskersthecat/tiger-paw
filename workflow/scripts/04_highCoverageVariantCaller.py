# -------------------------------
# 04_highCoverageVariantCaller.py
# tiger-paw
# Peter Reifenstein
# 
# Use: find variants in repeat segments
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

from collections import OrderedDict
import json

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("sam_file", help = "alignments to single reference")
parser.add_argument("reference_file", help = "reference sequence (.fa)")

parser.add_argument("--min_count", "-c", help = "report variants supported by this many reads", default = 2)
parser.add_argument("--min_depth", "-d", help = "minimum depth to report variant", default = 10)
parser.add_argument("--min_frequency", "-f", help = "minimum fraction to report variants, overwrites --min_count option", default = 0)

parser.add_argument("--REPEAT_LEN","-rlen", help = "reference represents repeated sequence with length rlen, adding two columns to output table", default = 0)

parser.add_argument("--MIN_COUNT_GIANT_INDEL","-mcg", help = "minimum count for giant indels to be automatically reported", default = 1)
parser.add_argument("--GIANT_INDEL_LEN","-lg", help = "minimum length of giant indels", default = 20)

parser.add_argument("--MODE","-m", help = "[SNP/INDEL] will only report variants on these types", default = "BOTH")
parser.add_argument("--KEEP_SUPPLEMENTARY","-s2", help = "don't skip supplementary alignments", action = "store_true")
parser.add_argument("--SKIP_SECONDARY","-s1", help = "skip secondary alignments", action = "store_true")


args = parser.parse_args()

sam_file = open(args.sam_file,'r')
ref_file = open(args.reference_file,'r')
min_count = int(args.min_count)
min_coverage = int(args.min_depth)
min_frequency = float(args.min_frequency)
variants_file = open(args.sam_file + ".VAR.tab",'w')
keep_supplementary = args.KEEP_SUPPLEMENTARY
skip_secondary = args.SKIP_SECONDARY

REPEAT_LEN = int(args.REPEAT_LEN)

header = "\"Reference Position\"	\"Type\"	\"Length\"	\"Reference\"	\"Allele\"	\"Count\"	\"Coverage\"	\"Frequency\""
if(args.REPEAT_LEN != 0):
   header += 	"	\"Segment_Number\"	\"Segment_Position\""

MIN_COUNT_GIANT_INDEL = args.MIN_COUNT_GIANT_INDEL
GIANT_INDEL_LEN = args.GIANT_INDEL_LEN

variants_file.write(header + "\n")

line = ref_file.readline()
ref_seq = ref_file.readline().strip()

ref_len = len(ref_seq)
print("[INFO] reference length:", ref_len)

depths = [0] * ref_len

Variations = OrderedDict({}) # Variations[variant_type][location][variation]

VARIANT_TYPES = ["SNV","INS","DEL"]

def incrementCount(variation_name, location, string):
    location = str(location)
    if string in Variations[variation_name][location].keys():
        Variations[variation_name][location][string] += 1
    else:
        Variations[variation_name][location][string] = 1

def formatPercentage(string):
    percent = str(string)
    nums = percent.split('.')
    beg = nums[0]
    while(len(beg) < 3):
        beg = "0" + beg 
    return beg + "." + nums[1][:1]

try:
    print("[INFO] Checking for json files on disc...")
    json_file = open(args.sam_file + ".VAR.tab.json", "r")
    Variations = json.load(json_file)
    json_file = open(args.sam_file + ".VAR.tab.depths.json", "r")
    depths = json.load(json_file)
    print("[INFO] Successfully loaded json file into memory")

except Exception as e:
    print("[INFO] No json file found, parsing variants ...")
    print("[INFO] Initializing memory ...")
    for variant_type in VARIANT_TYPES:
        Variations[variant_type] = OrderedDict({})

    for i in range (0, ref_len):
        for variant_type in VARIANT_TYPES:
            Variations[variant_type][str(i)] = {}
        
    print("[INFO] Reading alignments ...")

    total_aln = 0

    for line in sam_file.readlines():
        line = line.split()

        if(line[0][0] == '@'): # skip header fields:
            continue

        name = line[0]
        flag = line[1]

        if (not(keep_supplementary) and (int(flag) & 2048 == 2048)):
            print("[Find Variants] Skipping supplementary alignment with flag ", flag)
            continue

        if (skip_secondary and (int(flag) & 256 == 256)):
            print("[Find Variants] Skipping secondary alignment with flag ", flag)
            continue
        
        cigar = line[5]
        seq = line[9]
        start_pos = int(line[3]) - 1

        total_aln += 1
        if(total_aln % 100 == 0):
            print("Processed ", total_aln, " alignments...")

        # Generate the ALIGNED SEQUENCE
        align_seq = ""
        for i in range(0, int(start_pos)):   # if the alignment does not begin at position 0, buffer with N tokens:
            align_seq += "N"
        cur_num = ""
        pos = 0     # position in read sequence
        ref_pos = start_pos # position in reference sequence. should match with the length of align_seq

        for token in cigar:
            if token.isnumeric():
                cur_num += token
                continue
            num = int(cur_num)
            new_pos = pos
            if token == "M" or token == "=" or token == "X":
                new_pos = pos + num
                if token == "X":
                    incrementCount("SNV", ref_pos, seq[pos:new_pos])
                
                for i in range(num):
                    try:
                        depths[ref_pos + i] += 1
                    except:
                        print("Could not increment depth for position:", pos , "for alignment", line[:5])
                        exit()
                
                align_seq += seq[pos:new_pos]
                ref_pos = ref_pos + num

            elif token == "I":
                new_pos = pos + num
                incrementCount("INS", ref_pos, seq[pos:new_pos])
            
            elif token == "D":
                for i in range(0, num):
                    align_seq += "-"
                    depths[ref_pos + i] += 1
                
                incrementCount("DEL", ref_pos, ref_seq[ref_pos:num + ref_pos])
                ref_pos = ref_pos + num

            elif token == "S":
                new_pos = pos + num

            elif token == "H":
                pass
            else:
                print("CIGAR PARSE ERROR:: UNKNOWN TOKEN:", token)
                exit()

            pos = new_pos
            cur_num = ""

    print("[INFO] Dumping all variant info to file ...")

    with open(args.sam_file + ".VAR.tab.json", "w") as json_file:
        json.dump(Variations, json_file)

    with open(args.sam_file + ".VAR.tab.depths.json", "w") as json_file:
        json.dump(depths, json_file)

mode = args.MODE

min_len_INS_DEL = 1

if(min_frequency == 0):
    if(mode != "INDEL"):
        print("[CALL] Counting SNP and MNV which occur in at least ", min_count, " reads")
    if(mode != "SNP"):
        print("[CALL] Counting small (length  1) INS and DEL which occur in at least ", min_count, " reads")
        print("[CALL] Counting large (length >1) INS and DEL which occur in at least ", min_count, " reads")
        print("[CALL] Counting giant (length >", GIANT_INDEL_LEN ,") INS and DEL which occur in at least ", MIN_COUNT_GIANT_INDEL, " reads")

else:
    if(mode != "INDEL"):
        print("[CALL] Counting SNP and MNV which occur in at least ", min_frequency, " fraction of reads")
    if(mode != "SNP"):
        print("[CALL] Counting small (length  1) INS and DEL which occur in at least ", min_frequency, " fraction of reads")
        print("[CALL] Counting large (length >1) INS and DEL which occur in at least ", min_frequency, " fraction of reads")
        print("[CALL] Counting giant (length >", GIANT_INDEL_LEN ,") INS and DEL which occur in at least ", MIN_COUNT_GIANT_INDEL, " fraction of reads")

for location in range(len(depths)):
    
    depth = depths[int(location)]
    location = str(location)
    for variant_type in VARIANT_TYPES:
        if mode == "INDEL":
            if(variant_type == "SNV" or variant_type == "MNV"):
                continue
        
        if mode == "SNP":
            if(variant_type == "INS" or variant_type == "DEL"):
                continue

        if REPEAT_LEN != 0:
            repeat_number = int(int(location) / REPEAT_LEN) + 1
            repeat_pos = int(location) % REPEAT_LEN

        if location not in Variations[variant_type]:
            continue
        
        num_variations = 0
        for variation in Variations[variant_type][location]:
            alternate_count = Variations[variant_type][location][variation]

            if(min_frequency != 0):
                min_count = int(min_frequency * depth)

            varlen = len(variation)
            if (((varlen > GIANT_INDEL_LEN) and (alternate_count >= MIN_COUNT_GIANT_INDEL)) or (alternate_count >= min_count)):
                num_variations += 1
                freq = 0.0
                if(depth != 0):
                    freq = 100 * alternate_count / depth

                freq_string = formatPercentage(freq)
                varname = variant_type
                ref = ref_seq[int(location) : int(location) + varlen]
                alt = variation
                if varname == "SNV" and varlen > 1:
                    varname = "MNV"
                if varname == "INS":
                    ref = ""
                    for i in range(0, len(variation)):
                        ref += "-"
                    alt = variation
                if varname == "DEL":
                    ref = variation
                    alt = ""
                    for i in range(0, len(variation)):
                        alt += "-"
                
                # write this to the file
                variants_file.write(str(int(location) + 1) + "\t" + varname + "\t" + str(varlen) + "\t" + ref + "\t" + alt + "\t" + str(alternate_count) + "\t" + str(depth) + "\t" + freq_string)
                if(REPEAT_LEN != 0):
                    variants_file.write("\t" + str(repeat_number) + "\t" + str(repeat_pos))
                variants_file.write("\n")

        
        if variant_type == "DEL" and num_variations == 0:
            if (depth < min_coverage):
                variants_file.write(str(int(location) + 1) + "\t" + "LOW" + "\t" + "1" + "\t" + "-" + "\t" + "-" + "\t" + "0" + "\t" + str(depth) + "\t" + "000.0")
                if(REPEAT_LEN != 0):
                    variants_file.write("\t" + str(repeat_number) + "\t" + str(repeat_pos))
                variants_file.write("\n")

print("[INFO] Done")
