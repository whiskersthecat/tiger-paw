# -------------------------------
# 09_annotateSegments.py
# tiger-paw
# Peter Reifenstein
# 
# Use: label segment boundaries and assign variants to segment
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("gff_file", help = "sorted annotations to label (.gff)")
parser.add_argument("--ranges", "-r", help = "comma seperated list of tandem repeat cluster coordinates, e.g. Chromosome1:30000-31000,Chromosome8:4000-6000")
parser.add_argument("--shift", "-s", help = "shift all annotations by this amount", default = 0)
args = parser.parse_args()

gff_file = open(args.gff_file, 'r')
shift = int(args.shift)
output_gff_file = open(args.gff_file + ".labeled.gff", 'w')

def annotateSegment(line, segment, len):
    if(len != 0):
        line[8] = "LEN=" + str(len) + ";" + line[8] 
    if(segment != 0):
        line[8] = "SEGMENT=" + str(segment) +  ";" + line[8] 

    output_gff_file.write("\t".join(line) + "\n")

print("[info] Reading gff file")
current_segment = 0
previous_segment_coord = -1
for full_line in gff_file.readlines():
    if(full_line[0] == "#"):
        output_gff_file.write(full_line)
        continue
    
    full_line = full_line.strip()
    line = full_line.split('\t')
    chr = line[0]
    attributes = line[8].split(";")
    coords = (int(line[3]) + shift, int(line[4]) + shift)
    len = str( coords[1] - coords[0] )
    segment_id = -1
    new_segment = False
    new_attributes = attributes.copy()
    for i, attribute in enumerate(attributes):
        id, value = attribute.split("=")
        new_attributes[i] = attribute + ";"
        if id == "SEGMENT":
            segment_id = int(value)
        if id == "ORIGINAL_COORDINATES" or id == "ORIGINAL_COORD":
            new_attributes[i] = ""

    line[8] = "".join(new_attributes)
    if(segment_id == -1):
        annotateSegment(line, current_segment, 0)
    else:
        if(segment_id != current_segment):
            # new segment discovered
            if(previous_segment_coord != -1):
                annotateSegment([chr, "tiger-paw", "Segment", str(previous_segment_coord), str(coords[0]), ".", line[6], ".", ""], current_segment, str(coords[0] - previous_segment_coord))
            current_segment = segment_id
            previous_segment_coord = coords[0]

        annotateSegment(line, 0, 0)

annotateSegment([chr, "tiger-paw", "Segment", str(previous_segment_coord), str(coords[1]), ".", line[6], ".", ""], current_segment, str(coords[1] - previous_segment_coord))

print("[info] Done")

