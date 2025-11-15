# -------------------------------
# 10_segmentDrawer.py
# tiger-paw
# Peter Reifenstein
# 
# Use: draw variants in each segment, and the occurence of each segment type across the contig
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse
from PIL import Image, ImageDraw, ImageFont
import os

parser = argparse.ArgumentParser()
parser.add_argument("variants", help = "variants gff file, with segments and variants (from step 9 of tiger-paw)")
parser.add_argument("--repeat_cluster_sizes", "-rs", help = "provide list of repeat clusters and sizes used in naming, e.g. 'Chr1:10000,Chr8:8000' indicates segments beginning with Chr1 belong to one cluster with 10000 repeats", default = "")

args = parser.parse_args()

variants_file = open(args.variants, 'r')

variants_trace = open(args.variants + ".varianttraces", 'w')

MAX_SIZE = 10000

SCALE = 1

COL_BKGR = (255, 255, 255)

repeat_cluster_sizes = args.repeat_cluster_sizes.split(",")
repeat_cluster_size_dict = {}
try:
    for val in repeat_cluster_sizes:
        chr, size = val.split(":")
        repeat_cluster_size_dict[chr] = int(size)
except:
    print("[Error] repeat clusters sizes have not been provided correctly, see --help")

print("[draw] creating image object")


GIANT_MODEL_SCALE = 100 * SCALE

fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", int(20 * SCALE))
big_fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", int(80 * SCALE))
seq_fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", int(2 * SCALE))

model_fnt = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", int(1 * GIANT_MODEL_SCALE))


# visualization of segment CLADES
img = Image.new("RGB", (int(MAX_SIZE * 1.5 * SCALE), int(12000 * SCALE)), COL_BKGR)
draw = ImageDraw.Draw(img)

repeat_cluster_img_dict = {}
repeat_cluster_draw_dict = {}


COL_SEGMENT = (248, 244, 245) #(223, 230, 229) #(247, 228, 236)
COL_rRNA    = (229, 222, 224) #(203, 212, 211) #(230, 202, 205)
COL_SPACER  = (238, 232, 234) #(190, 200, 199)

COL_LABEL  = (66, 64, 65)

OUTLINE_SEG = (175, 175, 175)

REPRESENTATIVE_STROKE = 4

PRINT_LABELS = False


COL_SNV     = (0, 180, 161) #(150, 15, 103)
COL_MNV     = (0, 165, 180) #(150, 15, 103)

# COL_SNV     = (255, 0, 0)
COL_INS     = (255, 94, 151) #(224, 202, 230)

COL_RET     = (2, 105, 164) #(229, 255, 0)


# COL_INS     = (0, 255, 0)
COL_DEL     = (254, 124, 0) #(0, 0, 255)
COL_LARGE_DEL = (254, 216, 0) #(0, 0, 255)

COL_DEBUG   = (0, 0, 0)


COL_HIGHLIGHT = (255, 216, 50)
BLOAT = 1 * SCALE
HIGHLIGHT_BLOAT = 2 * SCALE

RETROTRANSPOSON_LENGTH = 1000

LARGE_DELETION_LENGTH = 500

OFFSET_COLORSHIFT = 1

RECT_HEIGHT = 15 * SCALE
SPACING_SCALE = 1.4
CORNER = (120 * SCALE, 40 * SCALE)

COLOR_RECT_LENGTH = 45 * SCALE
COLLAPSE_DELETIONS = True


for trc in repeat_cluster_size_dict:
   size = repeat_cluster_size_dict[trc]
   repeat_cluster_img_dict[trc] = Image.new("RGB", (int(RECT_HEIGHT * GIANT_MODEL_SCALE), int(size * GIANT_MODEL_SCALE)), COL_BKGR)
   repeat_cluster_draw_dict[trc] = ImageDraw.Draw(repeat_cluster_img_dict[trc])


n = 0

def convertColor(col, val):
   outputColor = tuple([int(a * val) for a in col])
   return outputColor

x, y = 0, 0
segment_coordinates = (0, 0)
prev_chrom = ""
compress = 1
deletion_shifts = []
seq = []
haplotype = ""

REPRESENTATIVE_SCALE = 7
HAPLOTYPE_CLUSTERING_SPACING = 55

rep_segment_n = -1
best_rep_count = 0
PRINT_REP = False
last_line = False

clade_stats = {}
cur_clade = ""

total_sum = 0
total_uniq = 0

# 1. Analysis

while True:

   PRINT_REP = False
   line = variants_file.readline()
   
   if not line:
      PRINT_REP = True
      last_line = True
   
   if(line[0:2] == "##"):
      clade_name = line[3:].strip()
      print("[clade]", clade_name)
      # draw the most common segment from before
      n += HAPLOTYPE_CLUSTERING_SPACING
      clade_stats[clade_name] = [0, 0]
      cur_clade = clade_name

      best_rep_count = 0
      if(rep_segment_n != -1):
         PRINT_REP = True
      # continue

   if(line[0:2] == "# " or PRINT_REP):

      if(PRINT_REP):
         rep = open("rep.txt", 'r')
         line = rep.readline()

      tmp = open("tmp.txt", "w")

      deletion_shifts = []

      n += 1

      if(rep_segment_n == -1):
         rep_segment_n = n

      fields = line[2:].split('\t')
      haplotype = fields[0]
      color = [int(col) for col in fields[1].split(",")]
      segments = [seg for seg in fields[2].split(" ")]

      variants_trace.write("\n>" + haplotype + "." + str(n) + "\nA")

      seq = list(fields[3])

      segment_location_dict = {trc : [] for trc in repeat_cluster_size_dict}

      new_best = False
      total_count = len(segments)
      if(total_count >= best_rep_count and not PRINT_REP):
         best_rep_count = total_count
         rep_write = open("rep.txt", "w")
         new_best = True
         rep_write.write(line)
         # print("new best segment: ", line)

      if(not PRINT_REP):
         clade_stats[cur_clade][0] += 1
         clade_stats[cur_clade][1] += total_count
         total_sum += total_count
         total_uniq += 1

      for seg in segments:
         if(seg[0] == ">"):
            seg = seg[1:]
         trc, location = seg.split("_")
         location = int(location[1:])
         try:
            segment_location_dict[trc].append(location)
         except:
            print("error: could not find tandem repeat cluster name", trc, " in list of tandem repeat clusters")
      

      ### DRAW the segment diagram giant model rectangles on the top
      if (not PRINT_REP):
         rx1, ry1 = 0, 0
         rx2 = rx1 + (RECT_HEIGHT * GIANT_MODEL_SCALE)

         for trc in repeat_cluster_size_dict:

            for seg_location in segment_location_dict[trc]:
               # seg_location = seg_location - 1
               sx1, sy1, sx2, sy2 = rx1 , ry1 + ((seg_location - 1) * GIANT_MODEL_SCALE), rx2 , ry1 + ((seg_location) * GIANT_MODEL_SCALE)

               repeat_cluster_draw_dict[trc].rectangle([sx1, sy1, sx2, sy2], fill = (color[0], color[1], color[2]))

               tx1, ty1 = sx1, sy1 + (GIANT_MODEL_SCALE * 0)

               repeat_cluster_draw_dict[trc].text((tx1, ty1), haplotype, font=model_fnt, fill = (255, 255, 255))

      ### DRAW the color label on the left side
      
      X, Y1, Y2 = CORNER[0], CORNER[1] + (RECT_HEIGHT * SPACING_SCALE * n), CORNER[1] + (RECT_HEIGHT * SPACING_SCALE * n) + RECT_HEIGHT
      if(PRINT_REP):
         X, Y1, Y2 = CORNER[0], CORNER[1] + (RECT_HEIGHT * SPACING_SCALE * (rep_segment_n - REPRESENTATIVE_SCALE)), CORNER[1] + (RECT_HEIGHT * SPACING_SCALE * (rep_segment_n - REPRESENTATIVE_SCALE)) + RECT_HEIGHT * SPACING_SCALE * (REPRESENTATIVE_SCALE - 1)

      x1, y1, x2, y2 = X, Y1, X + COLOR_RECT_LENGTH, Y2

      draw.rectangle([x1, y1, x2, y2], fill = (color[0], color[1], color[2]))
      draw.text((x1, y1), haplotype, font=fnt, fill = (255, 255, 255))


      ### DRAW the counts on the left side
      if(not PRINT_REP):
         x1, y1 = CORNER[0] * 0.6, Y1
         draw.text((x1, y1), str(total_count), font=fnt, fill = (color[0], color[1], color[2]))

      #### PREPROCESS AND FIND DELETION SHIFTS
      while(True):
         if(PRINT_REP):
            full_line = rep.readline()
            # print(full_line)
         else:
            full_line = variants_file.readline()
         
         if(full_line == "\n" or not full_line):
            break

         line = full_line.split('\t')
         try:
            start, end = int(line[3]), int(line[4])
            start, end = int(line[4]), int(line[3])
         except:
            print("error", line)
            exit()
         length = abs(end - start) + 1
         annotation_type = str(line[2])

         if(annotation_type == "Segment" or len(annotation_type) == 11):
            print("[segment]")
            deletion_shifts = length * [0]

            segment_coordinates = (start, end)

         if (annotation_type == "Variant"):
            position = (abs(start - segment_coordinates[0]))
            try:
               attributes = {str.split("=")[0]:str.split("=")[1] for str in line[8].strip().strip(";").split(";")}
            except Exception as e:
               print(e)
               print(line[8].strip().strip(";"))
               exit()
            vartype = attributes["TYPE"]
            
            ### SHIFT DELETIONS
            if(vartype == '"DEL"'):
               for i in range(position + 1, len(deletion_shifts)):
                  deletion_shifts[i] += length

            ### SHIFT LARGE INSERTIONS
            if(vartype == '"INS"' and length > RETROTRANSPOSON_LENGTH):
               for i in range(position + 1, len(deletion_shifts)):
                  deletion_shifts[i] -= length
            
         tmp.write(full_line)
         if(new_best):
            rep_write.write(full_line)
            # print("writing full line to rep:", full_line)

      tmp = open("tmp.txt", "r")
      #### PROCESS AND DRAW THE SEGMENT
      for line in tmp.readlines():

         line = line.split('\t')
         start, end = int(line[3]), int(line[4])
         start, end = int(line[4]), int(line[3])
         original_length = abs(end - start) + 1
         annotation_type = str(line[2])
         chrom = line[0]

         if(chrom != prev_chrom):
            prev_chrom = chrom
            print("[chrom]", chrom)
            n += 10
            

         start_pos = (abs(start - segment_coordinates[0]))
         end_pos = (abs(end - segment_coordinates[0]))

         try:
            if(deletion_shifts[start_pos] != 0):
               start_pos += deletion_shifts[start_pos]

            if(deletion_shifts[end_pos] != 0):
               end_pos += deletion_shifts[end_pos]
         
         except:
            print("[warning] failed to get deletion shift at position(s): ", start_pos, end_pos)

         length = abs(end_pos - start_pos) + 1

         color_offset = 1.0
         # if(n % 2 == 0):
         #    color_offset = OFFSET_COLORSHIFT\

         stroke_width = SCALE
         if(PRINT_REP):
            stroke_width = REPRESENTATIVE_STROKE * SCALE

         if(annotation_type == "Segment" or len(annotation_type) == 10):
            print("[segment] n = ", n, "length = ", length, "haplotype = ", haplotype)

            # length = abs(end - start) + 1 + deletion_shifts[end_pos]

            if(length > MAX_SIZE):
               compress = MAX_SIZE / length
               print("applying compression of ", compress)
            else:
               compress = 1

            X += COLOR_RECT_LENGTH * 2            
            
            x1, y1, x2, y2 = X + ((start_pos) * SCALE * compress) , Y1, X + ((end_pos) * SCALE * compress), Y2

            segment_coordinates = (start, end)
            draw.rectangle([x1, y1, x2, y2], fill = convertColor(COL_SEGMENT, color_offset), outline = OUTLINE_SEG, width = stroke_width)
         

         else:

            
            x1, y1, x2, y2 = X + ((start_pos) * SCALE * compress) , Y1, X + ((end_pos) * SCALE * compress), Y2
            if(PRINT_REP):
               x1 -= HIGHLIGHT_BLOAT
               x2 += HIGHLIGHT_BLOAT
            
            # col = (120, 245, 66)
            attributes = {str.split("=")[0]:str.split("=")[1] for str in line[8].strip().strip(";").split(";")}
            if("TYPE" in attributes):
               try:
                  attributes["TYPE"] = attributes["TYPE"].split('"')[1]
               except:
                  attributes["TYPE"] = attributes["TYPE"]

            if (annotation_type == "Variant"):

               x1 -= BLOAT
               x2 += BLOAT
               
               if("NEW" in attributes):
                  attributes["NEW"] = attributes["NEW"].split('"')[1]
               if("REF" in attributes):
                  attributes["REF"] = attributes["REF"].split('"')[1]
               # if("TYPE" in attributes):
               #    attributes["TYPE"] = attributes["TYPE"].split('"')[1]
               
               vartype = attributes["TYPE"]
               if(original_length > RETROTRANSPOSON_LENGTH and vartype == 'INS'):
                  # retrotransposon_shift += length - 1000
                  print("retrotransposon")
               
               outx1, outy1, outx2, outy2 = x1 - (HIGHLIGHT_BLOAT), y1 - (RECT_HEIGHT / 2), x2 + (HIGHLIGHT_BLOAT), y2 + (RECT_HEIGHT / 2)
               
               if(vartype == 'SNV'):
                  variants_trace.write(attributes["NEW"])
                  draw.rectangle([x1, y1, x2, y2], fill = COL_SNV)
               elif(vartype == 'MNV'):
                  variants_trace.write(attributes["NEW"])
                  draw.rectangle([x1, y1, x2, y2], fill = COL_MNV)

               
               elif(vartype == 'INS'):
                  variants_trace.write(attributes["NEW"])
                  if(original_length > RETROTRANSPOSON_LENGTH):
                     draw.rectangle([x1, y1, x1, y2], fill = COL_RET)
                     print("drawing text")
                     draw.text((x1 - (BLOAT), y1), "+" + str(original_length), font=fnt, fill = COL_RET)
                  else:
                     draw.rectangle([x1, y1, x2, y2], fill = COL_INS)
               
               elif(vartype == 'DEL'):
                  variants_trace.write(attributes["REF"])
                  start_pos = (abs(start - segment_coordinates[0]))
                  end_pos = (abs(end - segment_coordinates[0]))
                  
                  if(deletion_shifts[start_pos] != 0):
                     # print("shifted DELETION with deletion,", deletion_shifts[start_pos])
                     end_pos += deletion_shifts[start_pos]
                     start_pos += deletion_shifts[start_pos]
                  
                  # print("length of deletion:", end_pos - start_pos, start, end, start_pos, end_pos)
                  
                  x1, y1, x2, y2 = X + ((start_pos) * SCALE * compress) , Y1, X + ((end_pos) * SCALE * compress), Y2

                  x1 -= BLOAT
                  x2 += BLOAT
                  if(PRINT_REP):
                     x1 -= HIGHLIGHT_BLOAT
                     x2 += HIGHLIGHT_BLOAT

                  if(original_length > LARGE_DELETION_LENGTH):
                     draw.rectangle([x1, y1, x2, y2], fill = COL_LARGE_DEL)
                  else:
                     draw.rectangle([x1, y1, x2, y2], fill = COL_DEL)
                  
               else:
                  print("[warning] unknown variant type", vartype)
                  exit()
                  # draw.rectangle([x1, y1, x2, y2], fill = COL_DEBUG)
            

            if (annotation_type == "rRNA"):
               # continue
               try: 
                  draw.rectangle([x1, y1, x2, y2], fill = convertColor(COL_rRNA, color_offset), outline = OUTLINE_SEG, width = stroke_width)
               except Exception as e:
                  print(e)
                  print(x1, y1, x2, y2)
                  exit()
                  
               if(PRINT_REP and PRINT_LABELS):
                  draw.text((x1, y1 - (RECT_HEIGHT * 5)), attributes["TYPE"], font=big_fnt, fill = COL_LABEL)
               
            if (annotation_type == "spacer"):
               draw.rectangle([x1, y1, x2, y2], fill = convertColor(COL_SPACER, color_offset), outline = OUTLINE_SEG, width = stroke_width)
               if(PRINT_REP and PRINT_LABELS):
                  draw.text((x1, y1 - (RECT_HEIGHT * 5)), attributes["TYPE"], font=big_fnt, fill = COL_LABEL)
      
      tmp.close()

      if(PRINT_REP):
         rep_segment_n = n + 1
         rep.close()
      if(new_best):
         rep_write.close()

   if(last_line):
      break

# 2. Output
   
odirectory_name = "visualization"

if not os.path.exists(odirectory_name):
   print("[save] making output directory " + odirectory_name)
   os.mkdir(odirectory_name)

stats_file_name = args.variants + ".stats"
print("[save] writing stats to output file:", stats_file_name)

stats_file = open(stats_file_name ,'w')
for clade in clade_stats:
   count_uniq = clade_stats[clade][0]
   count_total = clade_stats[clade][1]
   percentage = 100 * count_total / total_sum
   stats_file.write("clade: " + str(clade) + "   unique: " + str(count_uniq) + "   total: " + str(count_total) + "    percentage: " + str(percentage) + "\n")

stats_file.write("total segments: " + str(total_sum) + "\n")
stats_file.write("total unique:   " + str(total_uniq) + "\n")

file_name = odirectory_name + "/output_clades.png" 
print ("[save] saving image " + file_name)
img.save(file_name, format = "png")

for trc in repeat_cluster_size_dict:
   file_name = odirectory_name + "/" + trc + ".png"
   print ("[save] saving image " + file_name)
   size = repeat_cluster_size_dict[trc]
   repeat_cluster_img_dict[trc].save(file_name, format = "png")
