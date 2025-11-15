# -------------------------------
# 04_phylogenyLabeler.py
# tiger-paw
# Peter Reifenstein
# 
# Use: label haplotypes based on grouping and divergence
# Reference: https://github.com/whiskersthecat/tiger-paw
# -------------------------------

import argparse
from collections import OrderedDict

import colorsys

parser = argparse.ArgumentParser()
parser.add_argument("tree_file", help = "phillips tree output")
parser.add_argument("--min_divergence", "-d", help = "minimum distance between nodes to create a new group (default 0.02)", default = 0.02)
parser.add_argument("--specific_divergence_list", "-l", help = "list of comma seperated divergence values to automatically add new group, e.g. 0.00435,0.21285,0.25308", default = "")
parser.add_argument("--specific_token_list", "-t", help = "list of comma seperated token strings to automatically add new group, e.g. jklmnor12345dehijqsuvy,jklmnor12345dehijqsuvyz", default = "")
parser.add_argument("--double_letter", "-z", help = "use two letters per haplotype ID", action = "store_true")

parser.add_argument("--add_colors", "-c", help = "assign colors to each haplotype", action = "store_true")
parser.add_argument("--custom_group_names", "-n", help = "names of each group, e.g. '14235'", default = "")
parser.add_argument("--custom_group_colors", "-r", help = "hues for each group (comma seperated), e.g. '221,33,2,95,144'", default = "")

args = parser.parse_args()

groups = "123456789abcdefghijklmnopqrstuvwxyz"
# groups = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnop"
cur_groups = ""
cur_group = 0
cur_group_name = ""

if (args.custom_group_names != ""):
    groups = args.custom_group_names


abridged_tokens = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
abridged_tokens = "abcdefghijklmnopqrstuvwxyz"

base = len(abridged_tokens)

group_counts = {}
group_haplotypes = OrderedDict({})

tokens = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

colors = ()

max_tokens = len(tokens)

cur_token = {}
cur_tokens = {}
cur_token[groups[0]] = 0
cur_tokens[groups[0]] = ""

INCLUDE_COLORS = args.add_colors

taxaPID = OrderedDict({})
taxaAbridgedName = OrderedDict({})

taxaColor = OrderedDict({})

prefix_size = 0
most_divergent_PID = ""
tree_file = open(args.tree_file, 'r')
names_file = open(args.tree_file + ".RENAMED.tsv", 'w')
stats_file = open(args.tree_file + ".STATS.txt", 'w')

if(len(cur_groups) == 0):
    # cur_groups += groups[0]
    # cur_group += 1
    # cur_group_name = groups[cur_group]
    # group_counts[cur_group_name] = 0
    # group_haplotypes[cur_group_name] = []


    cur_groups += groups[cur_group]
    cur_group_name = groups[cur_group]
    cur_token[cur_group_name] = 0
    group_counts[cur_group_name] = 0
    group_haplotypes[cur_group_name] = []

min_divergence = float(args.min_divergence)

specific_divergences = set()

if(args.specific_divergence_list != ""):
    for val in args.specific_divergence_list.strip().split(","):
        try:
            specific_divergences.add(float(val))
        except Exception as e:
            print(e)
            print(val, specific_divergences)
            print("[Error] divergence value not specified correctly")
            exit(0)

specific_tokens = set()

if(args.specific_token_list != ""):
    for val in args.specific_token_list.strip().split(","):
        try:
            specific_tokens.add((val))
        except Exception as e:
            print(e)
            print(val, specific_tokens)
            print("[Error] specific tokens not specified correctly")
            exit(0)

# -------------------------------
# 1. read tree file

# ignoreNextDivergence = False
            
NEW_GROUP = False

print("[INFO] reading tree file")

string = tree_file.read()
string = string[::-1]

mode = "none"

divergence = ""
name = ""

for char in string:
    if char == "\n" or char == ";":
        continue
    if mode == "number" :
        if char != ":" :
            divergence += char
            continue
        else:
            mode = "none"
            continue

    if mode == "name":
        if char == "," or char == "(":
            name = name[::-1]
            print("name:", name)
            if(len(name) > prefix_size):
                prefix_size = len(name)
            taxaPID[name] = cur_groups[len(cur_groups) - 1] + ":" + cur_tokens[cur_group_name]
            if(len(taxaPID[name]) > len(most_divergent_PID)):
                most_divergent_PID = taxaPID[name]

            if(name in specific_tokens):
                NEW_GROUP = True
            

            # if cur_group_name not in group_counts.keys():
            #     group_counts[cur_group_name] = 0 
            # if(INCLUDE_COLORS):
            #     taxaColor[name]= (curHue, curSaturation, curBrightness)

            # curSaturation -= 0.01
            # curBrightness -= 0.01
            
            
            num = group_counts[cur_group_name]
            # if num >= 27:
            #     num += 1
            secondLetter = abridged_tokens[ int(num / base) ]
            firstLetter = abridged_tokens [ num % base ]
            taxaAbridgedName[name] = cur_group_name + firstLetter
            if(args.double_letter):
                taxaAbridgedName[name] = cur_group_name + firstLetter + secondLetter
            group_counts[cur_group_name] += 1
            group_haplotypes[cur_group_name].append(name)

            mode = "none"
            name = ""
        
        else:
            name += char
            continue
    
    if char == "(":

        cur_tokens[cur_group_name] = cur_tokens[cur_group_name][:-1]
        cur_groups = cur_groups[:-1]

        print("intermediate tokens:", cur_tokens[cur_group_name])

        if (cur_tokens[cur_group_name] in specific_tokens):
            # print("Worked")
            # exit()
            NEW_GROUP = True
            specific_tokens.remove(cur_tokens[cur_group_name])

        # print("current_tokens:", cur_tokens[cur_group_name])
        if(len(cur_groups) > 0):
            cur_group_name = cur_groups[-1]
        else:
            cur_group_name = groups[0]
        continue

    
    if char == ")":
        ### FLUSH NUMBER as taxa divergence

        ### CHECK TAXA 
        if (divergence != ""):
            divergence = float(divergence[::-1])
        
            if(divergence > min_divergence or divergence in specific_divergences or cur_tokens[cur_group_name] in specific_tokens or NEW_GROUP == True):
                NEW_GROUP = False
                if (cur_tokens[cur_group_name] in specific_tokens):
                    specific_tokens.remove(cur_tokens[cur_group_name])
                
                print(" ## New group")
                if(group_counts[cur_group_name] > 0):
                    cur_group += 1
                    # curHue += 0.1
                    # curHue = curHue % 1
                    # curSaturation = 1
                    # curBrightness = 1
                
                cur_groups += groups[cur_group]
                cur_group_name = groups[cur_group]
                cur_token[cur_group_name] = 0
                group_counts[cur_group_name] = 0
                group_haplotypes[cur_group_name] = []

                print("### cur_group_name" + cur_group_name)
                
                cur_tokens[cur_group_name] = ""
            else:
                # maintain the current group
                cur_groups += cur_groups[len(cur_groups) - 1]

            print("Divergence for token " + cur_tokens[cur_group_name] + " = " + str(divergence))
            cur_tokens[cur_group_name] = cur_tokens[cur_group_name] + tokens[cur_token[cur_group_name] % max_tokens]
            # cur_tokens[cur_group_name] = cur_tokens[cur_group_name] + abridged_tokens[cur_token[cur_group_name]]

            cur_token[cur_group_name] += 1
            
            mode = "none"
            divergence = ""

        mode = "number"
        continue

    if char == ",":
        mode = "number"
        divergence = ""
        continue

    ### otherwise, start a new name:
    mode = "name"
    name = char


# -------------------------------
# 2. write the stats to a file

taxaColor = OrderedDict({})

curHue, curSaturation, curBrightness = 0, 0.8, 0.8

total_groups = len(group_haplotypes)
hue_adjust = 1 / total_groups


CUSTOM_COLORS = False
color_arr = ""
if(args.custom_group_colors != ""):
    try:
        color_arr = args.custom_group_colors.split(",")
        CUSTOM_COLORS = True
    except:
        print("Error: custom colors not specified correctly, see --help")

print("[INFO] writing stats")
# print(group_haplotypes)
n = 0
if(INCLUDE_COLORS):
    print("[INFO] assigning colors to each taxa")
    for groupName in group_haplotypes:
        if (CUSTOM_COLORS):
            try:
                curHue = float(color_arr[n])
            except:
                print("Warning: ran out of custom colors, using normal hue adjustment")
                CUSTOM_COLORS = False
        if (not CUSTOM_COLORS):
            curHue += hue_adjust
            curHue = curHue % 1
        curSaturation = 1
        curBrightness = 1
        count = group_counts[groupName]
        adjust = 0.4 / count
        for haplotype in group_haplotypes[groupName]:
            curSaturation -= adjust
            curBrightness -= adjust
            taxaColor[haplotype]= (curHue, curSaturation, curBrightness)
        
        n += 1

# print(taxaColor)

# print(colorsys.hsv_to_rgb(float(color_arr[0]), 1, 1))

stats_file.write("Appends Phylogenetic IDs (PIDs) for " + str(len(taxaPID.keys())) + " Taxa" + "\n")
stats_file.write("Most Divergent Branch: " + most_divergent_PID + "\n")

stats_file.write("Number of unique groups: " + str(len(group_counts.keys())) + "\n")
stats_file.write("Total Counts for each group: " + "\n")
for groupName in group_counts:
    stats_file.write("Group " + groupName + ": " + str(group_counts[groupName]) + "\n")


suffix_size = len(most_divergent_PID)

# -------------------------------
# 3. write the new names to a file

print("[INFO] writing new names")
# names_file.write("Taxa_New_Name\tOld_Name\n")
for name in taxaPID:
    suffix = taxaPID[name]
    while(len(suffix) < suffix_size):
        suffix += "-"
    prefix = name[14:17]
    new_name = prefix + "-" + suffix
    new_name = suffix
    
    names_file.write(name + "\t" + taxaAbridgedName[name] + "\t" +  new_name)

    if(INCLUDE_COLORS):
        rgbcolor = colorsys.hsv_to_rgb(taxaColor[name][0], taxaColor[name][1], taxaColor[name][2])
        names_file.write("\t" + str(int(255 * rgbcolor[0])) + "," + str(int(255 * rgbcolor[1])) + "," + str(int(255 * rgbcolor[2])))

    names_file.write("\n")




print("[INFO] done")
