# 4. Making Haplotypes

Overview:
By identifying high frequency variant sites in the reference segment, types of segments with specific variants, called a ***haplotype***, can be identified.

Using these, each segment in error-prone **long** reads can categorized into one of haplotypes.

In essence, this washes away every other base pair in the segment, allowing focus on true biological variation between repeats and ignoring almost all sequencing error in every read.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4a. Prepare folders and rename segments file.
```bash
mkdir alignments
mv reads/short.fa.chop.oriented.fa.segmented.fa reads/short.segments.fa
mv reads/long.fa.chop.oriented.fa.segmented.fa reads/long.segments.fa
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4b. Align all of the segments both datasets to the reference segment.

```bash
minimap2 -t *THREADS* --eqx -a reference/reference.fa reads/short.segments.fa > alignments/short.segments.sam
```
```bash
minimap2 -t *THREADS* --eqx -a reference/reference.fa reads/long.segments.fa > alignments/long.segments.sam
```
Note: the alignment for the long.segments won't be needed until part 5.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4c. Call variants on the alignment of the ***short.segments.sam*** alignment.

```bash
python3 04_highCoverageVariantCaller.py alignments/short.segments.sam reference/reference.fa --MODE SNP --min_frequency 0.02
```

Experiment with different min_frequency parameters for calling variants. \
It is recommended to set the min_frequency so that at least 10 biological segments in the tandem repeat cluster contain this variant.

You can estimate the total segments in the tandem repeat cluster:
```tsv
TOTAL_SEGMENTS_IN_TRC = total_segments_in_short_dataset / genomic_coverage_of_short_dataset
```
Then, min_frequency can be chosen accordingly:
```tsv
min_frequency = 10 / TOTAL_SEGMENTS_IN_TRC
```
For example, if TOTAL_SEGMENTS_IN_TRC = 1000 segments, this implies a --min_frequency of 0.01 requires at least 10 segments have this variant.


You may also call insertion and deletion variants. You may rename the output file to prevent overwriting.
```bash
python3 04_highCoverageVariantCaller.py alignments/short.segments.sam reference/reference.fa --MODE INDEL --min_frequency 0.02
```

Notes:
- Though insertions and deletion variants are common in repeat clusters, they are difficult to precisely measure in noisy long reads.
- Therefore, it is recommended to use simple SNV and MNV variants, with the option --MODE SNP.


It is also highly recommended to visualize variants to ensure that they make sense.\
If there are variants in homopolymer regions, they are probably not reliable to use (even if they may be real). You may opt to duplicate the output file and remove variants deemed this way.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4d. Find nucleotides at important locations in alignments.

Consider that the variants file from 4c contains these four variants:
```tsv
Reference Position	Type	Length	Reference	Allele	Count	Coverage	Frequency
500	SNV	1	C	T	678	30014	2.2589458252
1000	SNV	1	T	G	1305	30014	4.3479709468
1500	SNV	1	C	T	1313	30014	4.3746251749
2000	MNV	2	TC	GG	1301	30014	4.3346438328
```
For any given segment, the ***trace*** of the segment is the nucleotides at these 4 variant positions.

For example, based on the frequencies above, this is a good prediction of the most common trace among the 30014 segments:
```tsv
CTCTC
```
Notes:
- Although all these variant occur at very low frequency, importantly they are well within the boundaries of not being systematic sequencing errors. 
Therefore, they represent real, rare, biological polymorphism in this tandem repeat cluster.
- Any MNV (multi nucleotide variation) variants will take up multiple nucleotides in the trace, and are sorted to the end of the string.


Run this script.
```bash
python3 04_findVariants.py alignments/short.segments.sam alignments/short.segments.sam.VAR.tab reference/reference.fa
```
This will output a table which contains the nucleotide at each of the variant locations specified by the VAR file.

Next run this script.
```bash
bash 04_findTrace.sh alignments/short.segments.sam.variant_table.tab var1
```
Note:
- You can adjust **var1** to be another value for a different set of variants

This script will assemble the trace of each segment, and count the occurences of each trace.\

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4e. Creating haplotypes

```bash
mkdir haplotypes
```

Decide the minimum number of segments needed to support one trace for it to be real (MIN_SUPPORT) and considered as a haplotype.\
It is recommended to choose \*MIN_SUPPORT\* such that:
```tsv
MIN_SUPPORT >= genomic_coverage_of_short_dataset
```

```bash
awk '$1 >= *MIN_SUPPORT* {printf "HAPLOTYPE_%03dC%04d\t%s\n", NR, $1, $2 }' alignments/short.segments.sam.VAR.tab.Traces.Counts > haplotypes/short.segments.sam.VAR.tab.Traces.Counts.Named
```

Its very useful to rename the haplotypes based on their phylogenetic clustering.\
First convert to fasta format.
```bash
awk '{print ">" $1 "\n" $2}' alignments/short.segments.sam.VAR.tab.Traces.Counts.Named > alignments/short.segments.sam.VAR.tab.Traces.Counts.Named.fa
```
Now run this program and input the fasta file:
```tsv
https://www.ebi.ac.uk/jdispatcher/phylogeny/simple_phylogeny
```
Download the output as ".tree" format, and upload it into the alignments folder.\
Rename the file to "haplotypes.var1.tree"\
Run this script to map the original haplotype names to new haplotype names based on which class they cluster into.
```bash
python3 04_phylogenyLabeler.py alignments/haplotypes.var1.tree
```
This script will label each haplotype with an ID based their class. New classes will automatically be generated when there is a node with significant divergence. This threshold can be adjusted with --min_divergence to generate greater (or fewer) classes.



```bash
mkdir haplotypes
sort alignments/haplotypes.var1.tree.RENAMED.tsv > alignments/haplotypes.var1.tree.RENAMED.sorted.tsv
join alignments/haplotypes.var1.tree.RENAMED.sorted.tsv alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces.Counts.Named | awk '{print $2"\t"$4}' > haplotypes/haplotypes.var1.tsv
```

The haplotypes file should look something like this:
```tsv
1A	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGGGCCTCGAGATT
1D	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTAGCCTCGAGATT
1B	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGG---TCGAGATT
1A	CTCGCTAGTCAGAGTGAGTTCCTTAGATGACTTGGGG--TCGAGATT
7L	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTA---TCGAGATT
6D	CGTAACGAGACAGACAGTCTGTACGGCTAGTCCATGTCCGGTCTCAA
7I	CTCGCTAGTCAGAGTGAGTCCCTTAGATGACCTGTAG--TCGAGATT
...
```

Note:
- The "-" character indicates that the nucleotides is missing at the location
- Here, the first four haplotypes listed are all part of the same class (class 1)

Important:
- Haplotypes with many "-" characters can represent segments with large deletions at the beginning, middle, or end of the segment.
- If a haplotype has a single "-" character, it may not represent a real deletion, but rather a common misalignment due to sequencing errors. **In this case, the haplotype is not real, and is an imposter haplotype**. It is highly recommended to **manually remove these haplotypes from the haplotypes file**.
- If there are many haplotypes with many "-" characters at the *beginning or ending* of the segment, these could represent *truncated segments at the beginning or end of reads*. In this case, these are once again not real and imposter haplotypes. These haplotypes can be manually removed like above. **Alternatively, step 3c may be repeated with a more stringent (larger) value of -l**.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 4f. Creating another haplotype set

Repeat the entirity of step 4a-e, using a different variant frequency. Replace all occurences of **var1** with **var2**.

- It is highly recommended to make at least two distinct set of haplotypes, at different depths of variant frequency.
- For example, one at 0.01 frequency, and one at 0.02 frequency.
- Higher frequencies will produce tighter haplotypes with smaller traces, that are easier to categorize and more robust to sequencing errors.
- Lower frequencies will allow for rarer variations to be used, at the cost of decreased accuracy in haplotype assignment in part 5.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-05-BLOCKS](WORKFLOW-05-BLOCKS.md)