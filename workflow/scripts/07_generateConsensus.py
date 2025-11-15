# -------------------------------
# 07_generateConsensus.py
# tiger-paw
# Peter Reifenstein
# 
# Use: generate consensus sequence from variants file
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

from collections import OrderedDict

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("variantsfile", help = "table output from 04_highCoverageVariantCaller")
parser.add_argument("referencefile", help = "reference sequence (.fa)")

parser.add_argument("--annotations","-a", help = "comma seperated list of gff files to shift to match new sequence", default = "")
parser.add_argument("--exclude","-e", help = "comma seperated list of coordinates to ignore variants at", default = "")
parser.add_argument("--include","-i", help = "comma seperated list of coordinates to automatically include variants at", default = "")
parser.add_argument("--min_include_count","-mic", help = "minimum count to automatically include variants with include", default = 0)

parser.add_argument("--maxlength","-mlen", help = "cut any sequence and annotations after this point", default = 0)


parser.add_argument("--MIN_COUNT_ENORMOUS_INS","-mei", help = "minimum variant count to specially include enormous insertion", default = 10)
parser.add_argument("--MIN_COUNT_GIANT_INS","-mgi", help = "minimum variant count to specially include giant insertion", default = 6)

parser.add_argument("--MIN_BEST_COUNT_ENORMOUS_INS","-mbei", help = "minimum count of most common variant to specially include giant insertion", default = 1)
parser.add_argument("--MIN_BEST_COUNT_GIANT_INS","-mbgi", help = "minimum count of most common variant to specially include giant insertion", default = 3)


parser.add_argument("--MIN_COUNT_ENORMOUS_DEL","-med", help = "minimum variant count to specially include enormous deletion", default = 9)
parser.add_argument("--MIN_COUNT_GIANT_DEL","-mgd", help = "minimum variant count to specially include giant deletion", default = 6)

parser.add_argument("--MIN_BEST_COUNT_ENORMOUS_DEL","-mbed", help = "minimum count of most common variant to specially include enormous insertion", default = 1)
parser.add_argument("--MIN_BEST_COUNT_GIANT_DEL","-mbgd", help = "minimum count of most common variant to specially include giant insertion", default = 1)

parser.add_argument("--DONT_REMOVE_LOW_COVERAGE","-rl", help = "don't remove locations with LOW variation", action = "store_false")

parser.add_argument("--INCLUDE_BEST_INSERTION","-ibi", help = "if non-mutated is less than one half of reads, include best insertion (False/True)", action = "store_true")


args = parser.parse_args()

variants_file = open(args.variantsfile, 'r')
ref_file = open(args.referencefile, 'r')
log_file = open(args.variantsfile + ".consensus.log", 'w')
variants_used = open(args.variantsfile + ".consensus.variants", 'w')
lowcoverage_sites = open(args.variantsfile + ".consensus.lowcoverage", 'w')
new_seq_file = open(args.variantsfile + ".consensus.fa", 'w')

gff_annotation_file = open(args.variantsfile + ".consensus.gff", 'w')
gff_annotation_file.write("##gff-version 3\n")

remove_low_coverage = not(args.DONT_REMOVE_LOW_COVERAGE)
include_best_insertion = args.INCLUDE_BEST_INSERTION

min_include_count = int(args.min_include_count)

maxlength = int(args.maxlength)

name = ref_file.readline().strip()[1:]
ref_seq = ref_file.readline().strip()

# for mutated sequence
suffix = args.variantsfile.split('/')[-1]
new_seq_name = name + ".MUTATED." + suffix

ref_len = len(ref_seq)
log_file.write("[INFO] reference length:" + str(ref_len) + "\n")

read_names = ["Pos","Var","Len","Reference","Alternate","Frequency","Other_Count"]
Variations = {} #Variations[variant_type][location][variation] -> count

Variant_coordinates = {}

depths = [0] * ref_len
# Depth = {} #Depth[location] -> count
VARIANT_TYPES = ["SNV","MNV","INS","DEL","LOW"]
for variant_type in VARIANT_TYPES:
    Variations[variant_type] = OrderedDict({})

def log(msg):
    log_file.write(msg + "\n")
    log_file.flush()
    print(msg)

## 1. Read in variants file
log("[Mutate Reference] Reading in Variants")
line = variants_file.readline() # consume header
variants_used.write(line)

for full_line in variants_file.readlines():
    line = full_line.split('\t')
    var_type = line[1]
    location = int(line[0]) - 1
    variation = line[4]
    if(var_type == "DEL"):
        variation = line[3]

    count = int(line[5])
    depth = int(line[6])

    frequency = float(line[7])
    if not (location in Variations[var_type]):
        Variations[var_type][location] = {}
    
    Variant_coordinates[location] = True
    Variations[var_type][location][variation] = (count, full_line)
    try:
        depths[location] = depth
    except Exception as e:
        print(e)
        print(location)
        print(full_line)
        print(ref_len)

        exit()

for variant_type in VARIANT_TYPES:
    log("[INFO] " + str(len(Variations[variant_type])) + " total putative " + variant_type + " locations")

## 2. Mutate reference sequence
log("[Mutate Reference] Mutating Reference")

exclusion_coords = []
try:
    exclusion_coords = [(int(x) - 1) for x in args.exclude.split(",")]
    for x in exclusion_coords:
        log("excluding coordinate: " + str(x + 1))
except:
    log("[Mutate Reference] No coordinates excluded")

include_coords = []
try:
    include_coords = [(int(x) - 1) for x in args.include.split(",")]
    for x in include_coords:
        log("automatically including coordinate: " + str(x + 1))
except:
    log("[Mutate Reference] No coordinates automatically included")


new_seq = list(ref_seq.strip())

# GIANT
GIANT_INDEL_LEN = 20

MIN_COUNT_GIANT_DEL = int(args.MIN_COUNT_GIANT_DEL)
MIN_BEST_COUNT_GIANT_DEL = int(args.MIN_BEST_COUNT_GIANT_DEL)

MIN_COUNT_GIANT_INS = int(args.MIN_COUNT_GIANT_INS)
MIN_BEST_COUNT_GIANT_INS = int(args.MIN_BEST_COUNT_GIANT_INS)


# ENORMOUS
ENORMOUS_INDEL_LEN = 120

MIN_COUNT_ENORMOUS_DEL = int(args.MIN_COUNT_ENORMOUS_DEL)
MIN_BEST_COUNT_ENORMOUS_DEL = int(args.MIN_BEST_COUNT_ENORMOUS_DEL)

MIN_COUNT_ENORMOUS_INS = int(args.MIN_COUNT_ENORMOUS_INS)
MIN_BEST_COUNT_ENORMOUS_INS = int(args.MIN_BEST_COUNT_ENORMOUS_INS)


log("MIN_COUNT_ENORMOUS_INS " + str(MIN_COUNT_ENORMOUS_INS))
log("MIN_COUNT_GIANT_INS " + str(MIN_COUNT_ENORMOUS_INS))

ENORMOUS_DEL_MAX_LEN = 9000

coord_map = [0] * ref_len
new_coordinate = 0

def write_Annotation(type, start, end, ref, new, original_coord):
    gff_annotation_file.write(new_seq_name + "\tConsensus\tVariant\t" + str(start + 1) + "\t" + str(end + 1) + "\t.\t+\t.\tTYPE=\"" + type + "\";REF=\"" + ref + "\";NEW=\"" + new + "\";ORIGINAL_COORD=" + str(original_coord + 1) + "\n")

for location in range(ref_len):
    if(location % 1000000 == 0):
        print("Starting processing coordinate ", location)

    if(maxlength != 0 and location >= maxlength):
        new_seq[location] = ""
        continue

    if(location in Variant_coordinates):
        for variant_type in VARIANT_TYPES:
            if location not in Variations[variant_type]:
                continue

            if variant_type == "LOW":
                if(remove_low_coverage == "True"):
                    new_seq[location] = ""
                lowcoverage_sites.write(Variations[variant_type][location]["-"][1])
                continue

            if location in exclusion_coords:
                log("Skipping variant at coordinate: " + str(location))
                continue
            
            max_count = 0
            
            best_variation = ""
            depth = depths[location]
            non_mutated = depth
            large_indel_count = 0
            max_len = 0

            for variation in Variations[variant_type][location]:
                count = Variations[variant_type][location][variation][0]
                non_mutated -= count
                varlen = len(variation)

                if(count >= max_count or ((varlen >= GIANT_INDEL_LEN)) ):
                    
                    if(count >= max_count):
                        max_count = count
                        best_variation = variation

                    if(varlen > GIANT_INDEL_LEN):
                        # log("Found enormous indel " + Variations[variant_type][location][variation][1][0:50])
                        large_indel_count += count

                    if(varlen > max_len):
                        max_len = varlen
                        
                        
                        # log("Found giant indel with high count, count " + str(count) + " " + Variations[variant_type][location][variation][1][0:50])
                        # log("... new max count, variation is " + str(max_count) + " " + Variations[variant_type][location][best_variation][1][0:50])
            
            indel_type = "Small"
            if(max_len > ENORMOUS_INDEL_LEN):
                if(variant_type == "DEL" and large_indel_count >= MIN_COUNT_ENORMOUS_DEL and max_count >= MIN_BEST_COUNT_ENORMOUS_DEL):
                    indel_type = "Enormous"
                if(variant_type == "INS" and large_indel_count >= MIN_COUNT_ENORMOUS_INS and max_count >= MIN_BEST_COUNT_ENORMOUS_INS):
                    indel_type = "Enormous"
            
            elif (max_len > GIANT_INDEL_LEN):
                if(variant_type == "DEL" and large_indel_count >= MIN_COUNT_GIANT_DEL    and max_count >= MIN_BEST_COUNT_GIANT_DEL):
                    indel_type = "Giant"
                if(variant_type == "INS" and large_indel_count >= MIN_COUNT_GIANT_INS    and max_count >= MIN_BEST_COUNT_GIANT_INS):
                    indel_type = "Giant"
            

            if( not ( (location in include_coords) and (max_count >= min_include_count) ) ):
                # if there are more unmutated depth than mutated, ignore this variation 
                # exception: giant indels can have low depth
                if(indel_type == "Small"):
                    if(include_best_insertion == "True" and variant_type == "INS"):
                        if((non_mutated / depth) > 0.5):
                            best_variation = ""
                        
                    elif((non_mutated) > max_count):
                        best_variation = ""                

                # don't count variations in low coverage
                if(max_count <= 3 and large_indel_count <= 3):
                    best_variation = ""

                # enormous deletions cannot be too big
                if(variant_type == "DEL" and max_len > ENORMOUS_DEL_MAX_LEN):
                    best_variation = ""

            else:
                log("Automatically including variant at coordinate: " + str(location + 1))
                log("I:                                                         Automatically adding variant " + str(max_count) + " " + Variations[variant_type][location][best_variation][1][0:50])
                include_coords.remove(location)
                print(include_coords)

            if(best_variation != ""):

                if(indel_type == "Enormous"):
                    log("E: Adding enormous indel " + str(max_count) + " " + Variations[variant_type][location][best_variation][1][0:50])

                if(indel_type == "Giant"):
                    # log("location is at " + str(location))
                    log("G:                                                         Adding giant indel " + str(max_count) + " " + Variations[variant_type][location][best_variation][1][0:50])

                variants_used.write(Variations[variant_type][location][best_variation][1])

                if variant_type == "SNV":
                    new_seq[location] = best_variation

                    write_Annotation("SNV", new_coordinate, new_coordinate, ref_seq[location], best_variation, location)
                    # gff_annotation_file.write("Ls_NOR_Ref\tConsensus\tVariant\t1\t1574	.	+	.	TYPE="SNV",REF="",NEW="",ORIGINAL_COORD=location\n")

                if variant_type == "MNV":
                    for i, nucleotide in enumerate(best_variation):
                        new_seq[location + i] = nucleotide
                    write_Annotation("MNV", new_coordinate, new_coordinate + (len(best_variation) - 1), ref_seq[location : location + len(best_variation)], best_variation, location)
                
                if variant_type == "INS":
                    new_seq[location] = best_variation + new_seq[location]
                    write_Annotation("INS", new_coordinate, new_coordinate + (len(best_variation) - 1), "", best_variation, location)

                if variant_type == "DEL":
                    # Delete this many bases
                    for i, nucleotide in enumerate(best_variation):
                        new_seq[location + i] = ""
                    write_Annotation("DEL", new_coordinate, new_coordinate + (len(best_variation) - 1), best_variation, "", location)

    coord_map[location] = new_coordinate
    # increment coordinate by length of new sequence
    try:
        new_coordinate += len(new_seq[location])
    except Exception as e:
        print(e)
        print("location:", location, " new_seq length", len(new_seq))
        exit()



new_seq_str = "".join(new_seq)
new_seq_file.write(">" + new_seq_name + "\n")
new_seq_file.write(new_seq_str + "\n")

## 3. Convert coordinates

gff_file_names = args.annotations.split(",")

for gff_file_name in gff_file_names:
    if(gff_file_name == ""):
        continue

    log("[Mutate Reference] Converting Coordinates in gff file " + gff_file_name)
    gff_file = open(gff_file_name, 'r')
    gff_converted_file = open(gff_file_name + "." + suffix + ".consensus.gff", 'w')

    # read through gff file and change coordinates
    # header = gff_file.readline()
    # gff_converted_file.write(header)

    header = False
    for full_line in gff_file.readlines():
        if(full_line[0] == "#"):
            if (header == False):
                gff_converted_file.write(full_line)
                header = True
            continue

        full_line = full_line.strip()
        line = full_line.split('\t')
        if(line[0] != name):
            continue

        coords = (int(line[3]), int(line[4]))
        try:
            new_coords = (str(coord_map[coords[0] - 1] + 1), str(coord_map[coords[1] - 1] + 1))
        except Exception as e:
            print(e)
            print(full_line)

        # check for deletion of  feature
        # for features with no length, check the following coordinates new location
        if(coords[0] == coords[1]):
            if(coord_map[coords[0] - 1] == coord_map[coords[1]]):
                log("[Deleted feature] " + full_line)
                continue

        elif(new_coords[0] == new_coords[1] or int(new_coords[1]) > len(new_seq)):
            log("[Deleted feature] " + full_line)
            continue

        line[0] = new_seq_name
        line[3] = new_coords[0]
        line[4] = new_coords[1]
        line[8] = line[8] + ";ORIGINAL_COORDINATES=(" + str(coords[0]) + "-" + str(coords[1]) + ")"

        gff_converted_file.write("\t".join(line) + "\n")

print("[INFO] Done")