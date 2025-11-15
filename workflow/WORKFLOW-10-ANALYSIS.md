# 10. Analysis

Note: This optional step contains tools used to generate graphics depicting the results of the Assembly of the Nucleolus Organizer Regions in Lettuce. However, it may require some modifications for other projects.

### 10a. K-mer completeness.

The incorporation of a correctly assembled TRC generates novel sequence.\
Goal: calculate the percent of kmers in any number of datasets (reads) that occur in databases (genomes).\

Prepare the following files:
- genome/genome.fa *Final genome.*
- genome/noTRC.fa. *Original genome without TRCs.*
- genome/Contig1.fa *Assembled TRCs*
- genome/originalContig1.fa *Original portion of noTRC.fa that represents TRC*

Install required software KMC.
```bash
mamba install kmc -c bioconda
```
Run script.
```bash
bash scripts/10_computeKmers.sh
```

### 10b. Exploring diversity of segments

> As discovered in part 4, tandem repeat clusters have many mutations distinguishing indiviual segments. Visualizing where these mutations occur, both along the TRC, and within each segment, allows for understanding of the evolution of the TRC.

Extract sequences of segments and assign names to each segment.
```bash
grep -P "\tSegment\t" genome/genome.tigerpaw.gff > genome/genome.tigerpaw.segments.gff
awk 'BEGIN {FS=OFS="\t"} {count[$1] +=1; $3 = sprintf("%s_s%05d", $1, count[$1]); print}' genome/genome.tigerpaw.segments.gff > genome/genome.tigerpaw.segments.named.gff
bedtools getfasta  -nameOnly -fi genome/genome.fa -bed genome/genome.tigerpaw.segments.named.gff -fo genome/genome.segments.fa
```

Split the gff file into one for each segment.
```bash
awk '$3 != "Segment"' genome/genome.tigerpaw.gff > genome/genome.tigerpaw.notsegments.gff
cat genome/genome.tigerpaw.notsegments.gff genome/genome.tigerpaw.segments.named.gff | sort -k4,4n -k2,2r > genome/genome.tigerpaw.named.gff
python3 10_segmentGFFWriter.py genome/genome.tigerpaw.named.gff
```

Convert to tab-delimited, and cluster identical segments.
```bash
awk '{printf("%s", $1); getline; printf("\t%s\n", $1);}' genome/genome.segments.fa > genome/genome.segments.tab
awk 'BEGIN {FS=OFS="\t"} {a[$2] = a[$2]$1" "; c[$2] = c[$2] += 1} END{for (i in a) print c[i] "\t" a[i] "\t" i}' genome/genome.segments.tab | sort -k1,1nr -k2,2 > genome/genome.segments.grouped.tab
```
Inspect the output file genome/genome.segments.grouped.tab. Each line lists the number of occurences of this segment sequence, which segments it occurs in, and the sequence itself.\
Assign preliminary names to each segment type.
```bash
awk 'BEGIN {FS=OFS="\t"} {printf("%s.%s\t%s\t%s\n", substr($2, 2, 11), $1, $2, $3)}' genome/genome.segments.grouped.tab > genome/genome.segments.grouped.named.tab
awk 'BEGIN {FS=OFS="\t"} {printf(">%s\n%s\n", $1, $3)}' genome/genome.segments.grouped.named.tab > genome/genome.segments.grouped.named.fa
```

Now cluster these unique segments.
```bash
clustalo --threads 20 --force -i genome/genome.segments.grouped.named.fa -o genome/genome.segments.grouped.named.clustalio.fa
```
```bash
clustalw -infile=genome/genome.segments.grouped.named.clustalio.fa -TREE -OUTPUTTREE=phylip
```

Assign labels to each segment.
```bash
python3 scripts/04_phylogenyLabeler2.py genome/genome.segments.grouped.named.clustalio.ph -d 0.02 -c --double_letter
```
Notes: 
- The -c option will assign a color to each segment type. Segments that cluster together will be given the same hue, but with varying saturation. This is useful in visualization.
- Adjust -d to allow for more clustering, or use --specific_token_list to force clustering at certain nodes

Relabel original segments.
```bash
awk 'BEGIN {FS=OFS="\t"} NR == FNR {map[$1] = $2"\t"$3; next;} { key = $1; if (key in map) { $0 = $0 "\t" map[key]; } print;}' genome/genome.segments.grouped.named.tab genome/genome.segments.grouped.named.clustalio.ph.RENAMED.tsv > genome/genome.segments.grouped.named.renamed.tab
```

Generate a GFF file with the relabeled segments
```bash
python3 scripts/10_segmentGFFMerger.py genome/genome.segments.grouped.named.renamed.tab
```
Notes:
- This file contains the exact variants of each unique segment type, the sequence of each segment, and the cluster it belongs to.

From this file, a visualization can be generated. First, install Python Imaging Library (Pillow).
```bash
mamba install pillow
```

Generate a visualization. Replace Chr1:50 with a comma seperated list of cluster sizes (in segments)
```bash
python3 scripts/10_segmentDrawer.py genome/genome.segments.grouped.named.renamed.tab.merged.gff -rs Chr1:50
```
