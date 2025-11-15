# -------------------------------
# 09_concatenateFasta.py
# tiger-paw
# Peter Reifenstein
# 
# Use: concatenate two fasta files while updating gff coordinates
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("seq1", help = "first sequence (.fa)")
parser.add_argument("seq2", help = "second sequence (.fa)")
parser.add_argument("--annotations","-a", help = "comma seperated list of gff files to shift to match new sequence")
parser.add_argument("--dontupdatecoords","-d", help = "do not update coordinates in provided gff files (only change reference names)", action = "store_true")
parser.add_argument("--customprefix","-p", help = "custom name for input sequence", default = "")
args = parser.parse_args()

update = not(args.dontupdatecoords)
seq1_file = open(args.seq1, 'r')
seq2_file = open(args.seq2, 'r')

fname1 = args.seq1.split('/')[-1]
fname2 = args.seq2.split('/')[-1]

name1 = seq1_file.readline()[1:].strip()
name2 = seq2_file.readline()[1:].strip()

seq1 = seq1_file.readline().strip()
seq2 = seq2_file.readline().strip()

seqname = fname1 + ".cat." + fname2

filename = seqname
customprefix = args.customprefix
if(customprefix != ""):
    output_seq_file = open(args.seq1 + ".cat." + fname2, 'w')
else:
    output_seq_file = open(args.seq1 + ".cat." + fname2, 'w')

output_seq_file.write(">" + seqname + "\n" + seq1 + seq2 + "\n")

shift = len(seq1)

gff_file_names = args.annotations.split(",")

for gff_file_name in gff_file_names:

    print("Converting Coordinates in gff file " + gff_file_name)
    gff_file = open(gff_file_name, 'r')
    new_gff_file_name = gff_file_name
    # print("Updating coordinates ")
    # if(update == True and customprefix == ""):
    #     new_gff_file_name = gff_file_name + ".cat." + name1 + ".gff"
    # else:
    new_gff_file_name = gff_file_name + ".cat.gff"
    gff_converted_file = open(new_gff_file_name, 'w')

    for full_line in gff_file.readlines():
        if(full_line[0] == "#"):
            gff_converted_file.write(full_line)
            continue
        
        full_line = full_line.strip()
        line = full_line.split('\t')
        coords = (int(line[3]), int(line[4]))
        new_coords = (str(coords[0] + shift), str(coords[1] + shift))

        line[0] = seqname
        if(update == True):
            line[3] = new_coords[0]
            line[4] = new_coords[1]

        gff_converted_file.write("\t".join(line) + "\n")
