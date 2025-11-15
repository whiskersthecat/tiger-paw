# -------------------------------
# 10_segmentGFFWriter.sh
# tiger-paw
# Peter Reifenstein
# 
# Use: organize GFF input file into one per segment
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("variants", help = "variants gff file, with segments and variants (from step 9 of tiger-paw)")
args = parser.parse_args()

variants_file = open(args.variants, 'r')

if not os.path.exists("segment_gff"):
   os.mkdir("segment_gff")


ofile  = open("segment_gff/" + "unplaced" + ".gff","w")

for full_line in variants_file.readlines():
   if(full_line[0] == "#"):
      continue
   line = full_line.split('\t')
   name = line[2]

   if(len(name) == 11):
      # indicates this is a segment
      ofile = open("segment_gff/" + name + ".gff","w")
   ofile.write(full_line)
