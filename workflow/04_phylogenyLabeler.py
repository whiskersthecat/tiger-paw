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

parser = argparse.ArgumentParser()
parser.add_argument("tree_file", help = "file output from https://www.ebi.ac.uk/jdispatcher/phylogeny/simple_phylogeny")
parser.add_argument("--min_divergence", "-d", help = "minimum distance between nodes to create a new group", default = 0.02)
parser.add_argument("--specific_divergence_list", "-l", help = "list of comma seperated divergence values to automatically add new group, e.g. 0.00435,0.21285,0.25308", default = "")

args = parser.parse_args()

groups = "123456789"
cur_groups = ""
cur_group = 0
cur_group_name = ""

abridged_tokens = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

base = len(abridged_tokens)

group_counts = {}

tokens = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()¡™£¢∞§¶§¶•ªºœ∑´®†¥¨¨ˆøπåß∂ƒ©˙∆∆˚¬≈ç√∫˜µ≤"

cur_token = {}
cur_tokens = {}
cur_token[groups[0]] = 0
cur_tokens[groups[0]] = ""

taxaPID = OrderedDict({})
taxaAbridgedName = OrderedDict({})

prefix_size = 0
most_divergent_PID = ""
tree_file = open(args.tree_file, 'r')
names_file = open(args.tree_file + ".RENAMED.tsv", 'w')
stats_file = open(args.tree_file + ".STATS.txt", 'w')

if(len(cur_groups) == 0):
    cur_groups += groups[0]
    cur_group += 1
    cur_group_name = groups[0]

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

# -------------------------------
# 1. read tree file

print("[INFO] reading tree file")
for line in reversed(tree_file.readlines()):
    # if the line is a :, add a new branch
    if(line[0] == ":"):
        # check the divergence of this branch to check for new group
        divergence = float(line[1:8])
        
        if(divergence > min_divergence or divergence in specific_divergences):
            print(" ## New group")
            cur_groups += groups[cur_group]
            cur_group_name = groups[cur_group]
            cur_token[cur_group_name] = 0
            cur_group += 1
            cur_tokens[cur_group_name] = ""
        else:
            # maintain the current group
            cur_groups += cur_groups[len(cur_groups) - 1]

        print("Divergence for token " + tokens[cur_token[cur_group_name]] + " = " + str(divergence))
        cur_tokens[cur_group_name] = cur_tokens[cur_group_name] + tokens[cur_token[cur_group_name]]
        cur_token[cur_group_name] += 1

    # if the line ends with a ), remove the last token added
    if(line[0] == "("):
        cur_tokens[cur_group_name] = cur_tokens[cur_group_name][:-1]
        cur_groups = cur_groups[:-1]
        # print("current_tokens:", cur_tokens[cur_group_name])
        if(len(cur_groups) > 0):
            cur_group_name = cur_groups[-1]
        else:
            cur_group_name = groups[0]
        continue

    # if the line contains the id of a taxa, which occurs when it does not start with : or ( add this taxa name to the dictionary
    if(line[0] != ":"):
        name = line.split(":")[0]
        print(name)
        if(len(name) > prefix_size):
            prefix_size = len(name)
        taxaPID[name] = cur_groups[len(cur_groups) - 1] + ":" + cur_tokens[cur_group_name]
        if(len(taxaPID[name]) > len(most_divergent_PID)):
            most_divergent_PID = taxaPID[name]
        if cur_group_name not in group_counts.keys():
            group_counts[cur_group_name] = 0 
        
        num = group_counts[cur_group_name]
        if num >= 27:
            num += 1
        firstLetter = abridged_tokens[ int(num / base) ]
        secondLetter = abridged_tokens [ num % base ]
        taxaAbridgedName[name] = cur_group_name + firstLetter + secondLetter
        taxaAbridgedName[name] = cur_group_name + secondLetter
        group_counts[cur_group_name] += 1

# -------------------------------
# 2. write the stats to a file

print("[INFO] writing stats")
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
    
    names_file.write(name + "\t" + taxaAbridgedName[name] + "\t" +  new_name + "\n")

print("[INFO] done")
