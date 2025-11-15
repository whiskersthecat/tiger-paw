# -------------------------------
# 10_segmentGFFMerger.sh
# tiger-paw
# Peter Reifenstein
# 
# Use: combine segment gff files into one
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("haplotypes", help = "renamed haplotype file from 04_phylogenyLabeler, with haplotype of each segment. First 10 characters of field one are looked up as file in folder segment_gff (from script spreadSegmentToGFF)")
args = parser.parse_args()

haplotype_file = open(args.haplotypes, 'r')

if not os.path.exists("segment_gff"):
   print("Error: run the script 10_segmentGFFWriter.py")
   exit()

ofile  = open(args.haplotypes + ".merged.gff", "w")

haplotype_group =""
prev_group = ""
for full_line in haplotype_file.readlines():
   line = full_line.split('\t')
   name = line[0]
   name = line[0][0:11]
   haplotype = line[1]
   haplotype_group = haplotype[0]
   try:
      color = line[3].strip()
   except:
      print("Error: add the color column (-c) when running 04_phylogenyLabeler")
      exit()

   try:
      segments = line[4].strip()
   except:
      print("Warning: no segments column provided")
      segments = ""

   try:
      seq = line[5].strip()
   except:
      print("Warning: no sequences column provided")
      seq = ""

   if(haplotype_group != prev_group):
      ofile.write("## " + haplotype_group + "\n")
   prev_group = haplotype_group

   path = "segment_gff/" + name[0:13] + ".gff"
   try:
      gff_file = open(path, 'r')
   except:
      print("File does not exist:" , path)
      print("Make sure to run 10_segmentGFFWriter.py first")
      exit()
   
   ofile.write("# " + haplotype + "\t" + color + "\t" + segments + "\t" + seq + "\n")
   for line in gff_file.readlines():
      ofile.write(line)
   ofile.write("\n")
