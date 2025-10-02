# 0. Introduction


These **definitions** are used throughout the workflow.

* ***segment***: one repeat copy of a tandem repeat. A single long read should span many segments.
* ***chimera***: a single read in which two different (unrelated) molecules are concatenated together.
* ***duplex***: a single read in which the same molecule is read twice, with forward and reverse strand concatenated together.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>




**tiger-paw** requires two large DNA sequencing datasets, ideally both with coverage at least 50x. For *long.fa*, the more the better.

* ***short.fa***. Reads are up to 30,000 base pairs and each spans at least one segment, and are **very accurate**. They are assumed to contain no chimeras or duplexes. For example, *Pacbio HiFi Reads*. 
> These reads will be used to generate common Haplotypes that will be used in assembly.


* ***long.fa***. Reads are longer than 30,000 base pairs and can span at **least three segments, ideally up to 20 or more**. They are assumed to be innacurate, contain duplexes, and may contain chimeras. For example, *Oxford Nanopore Reads*. 
> These reads will be used to actually assemble the tandem repeat cluster. 

Note: in the case there are multiple datasets of either type, it is best to keep them seperate and run them in parallel. This helps confirm phenomena that is replicated among all datasets.


<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Main Dependencies:
- minimap2 (https://github.com/lh3/minimap2)
- blast Version 2.16.0 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)
   - Note: different blast Versions have very different performance. Newest blast version 2.17.0 struggled in identifying duplexes in HiFi reads.
- bedtools
- Python3

You can install the dependencies using mamba:
```bash
mamba install -c bioconda minimap2 blast=2.16.0 bedtools python=3.11
```

Supplementary Software:\
It is highly recommended to use one of the following visualization software to visualize alignments on multiple steps
- Interactive Genomics Viewer (https://igv.org/)
- CLC Genomics Workbench (https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/qiagen-clc-genomics/?cmpid=QDI_GA_DISC_CLC_SA&gad_source=1&gad_campaignid=6942460474&gbraid=0AAAAADbyWl0yf7q1l0zt_ZecsD7avy_bK&gclid=Cj0KCQjwovPGBhDxARIsAFhgkwSB0_xW4BRg1U7Z2qqdnT1Qtfd1Hrt5fg4xcGbaECuH1LFYdof3tpEaAsMCEALw_wcB)

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 0A. Copy reads

Note: in all commands presented in the workflow, text in \* \* should be replaced with user parameters.

```bash
mkdir reads
cp *path_to_short_reads* reads/short.fa
cp *path_to_long_reads* reads/long.fa
```


# 1. Reference

```bash
mkdir reference
```


# 2. Reads

Overview: for each dataset, generate a highly refined subset of reads which contain repeat segments.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 2a. Find reads with BLAST hits to reference segment.

```bash
mkdir blastN
mkdir blastN/DB
```
```bash
makeblastdb -in reads/short.fa -out blastN/DB/short.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
makeblastdb -in reads/long.fa -out blastN/DB/long.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
```
Modify \*THREADS\* based on available computing resources.
```bash
blastn -task blastn -query reference/reference.fa -db blastN/DB/short.BDB -out blastN/reference.vs.short.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
blastn -task blastn -query reference/reference.fa -db blastN/DB/long.BDB -out blastN/reference.vs.long.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 2b. Reformat BLAST hits.

Collect large blast hits and reformat table.

For ***short*** dataset, the value \*SIZE\* should be adjusted to be roughly 60% of the reference segment length.
For ***long*** dataset, the value \*SIZE\* should be adjusted to be roughly 20% of the reference segment length.\

```bash
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.short.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.short.sorted.blastN
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.long.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.long.sorted.blastN
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 2c. Orienting of reads and finding duplexes.

In tandem repeats, segments may switch orientation at certain loci.
Notes:

- If this phenomenon is not observed in ***short*** dataset, this phenomenon does not exist. 
In this case, any switches in orientation observed in the ***long*** dataset represent duplex reads.
These duplex reads should be chopped, and the full length read should be removed from the dataset.

- On the other hand, if such a phenomenon is observed in the ***short*** dataset, these likely truly exist. There are two cases:
   + The orientation switch occurs in a single or few segments. In this case, duplex reads from the ***long*** dataset should be chopped like above.
   + The orientation switch occurs systematically at a certain breakpoint. In this case, distinguishing duplex reads in the long.fa dataset from real reads spanning this region is not possible at this point. Reads marked as duplex in ***long*** dataset should be used with caution during assembly.

Run the following script on both datasets.
```bash
python3 02_chopDuplexAndOrientReads.py reads/short.fa blastN/reference.vs.short.sorted.blastN --include_original
python3 02_chopDuplexAndOrientReads.py reads/long.fa blastN/reference.vs.long.sorted.blastN --include_original
```

This script will generate multiple output files.
1. reads/short.fa.**report** a count of how many reads were in each of 5 categories
2. reads/short.fa.**readids** the categorization of each read
3. reads/short.fa.**chop.oriented.fa** an output fasta file that includes modified reads:
   - With --include_original, any reads marked as duplexes will have two corresponding reads in the output file:
      - read_name: the modified read.
      - read_name_dp (duplex perfect) or read_name_dbf (duplex bias forward) or read_name_dbr (duplex bias reverse). This read has the sequence of the original read, and indicates its duplex categorization.
   - Otherwise, only modified reads are output

# 3. Segments

Overview: for each dataset, cut reads containing tandem repeats into segments.\
Consider this simplified example with a very short 5 base pair reference segment of 'AAATG':

**input_reads.fa:**
```tsv
read_01  G AAATG AAATG AATTG AAATG AAATG AA
read_02  CCCCGCCC AAATG A
```
**blastHits.blastN:**
```tsv
read_01  1 2---- 3---- 4---- 5---- 6---- 7-
read_02           1----
```
**output_reads.fa:**
```tsv
read_01_segment_02 AAATG
read_01_segment_03 AAATG
read_01_segment_04 AATTG
read_01_segment_05 AAATG
read_01_segment_06 AAATG
read_01_segment_07 AA

read_02_segment_01 AAATG
```

Notes:
- the first segment of read_01 was removed because it was small (only 1 base pair)
- read_02 might contain the linkage between the chromosome and the tandem repeat. Alternatively, it could be a chimeric read.
- read_01_segment_04 contains a variant (T instead of A). This read could be a variant and serve useful in assembly. Alternatively, it could be a sequencing error that will confound downstream assembly.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 3a. Redo blast hits on new dataset with oriented and chopped reads.


```bash
makeblastdb -in reads/short.fa.chop.oriented.fa -out blastN/DB/short.fa.chop.oriented.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
makeblastdb -in reads/long.fa.chop.oriented.fa -out blastN/DB/long.fa.chop.oriented.BDB -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids  
```
Modify \*THREADS\* based on available computing resources.
```bash
blastn -task blastn -query reference/reference.fa -db blastN/DB/short.fa.chop.oriented.BDB -out blastN/reference.vs.short.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
blastn -task blastn -query reference/reference.fa -db blastN/DB/long.fa.chop.oriented.BDB -out blastN/reference.vs.long.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 3b. Reformat BLAST hits.

Like step 2b, collect large blast hits and reformat table.

For both datasets, dataset, the value \*SIZE\* should be adjusted to be roughly 20% of the reference segment length.\
```bash
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.short.fa.chop.oriented.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.short.fa.chop.oriented.sorted.blastN
awk -F'\t' '$6 >= *SIZE*' blastN/reference.vs.long.fa.chop.oriented.blastN | sort -k17nr -k2,2 -k13,13n > blastN/reference.vs.long.fa.chop.oriented.sorted.blastN
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 3c. Segment the reads based on the BLAST hits.

Run the following script on both datasets.\

The value for \*SIZE\* determines the tolerance to include the segments at the edges of the reads.\
A larger (stricter) value will produce fewer segments.\
A smaller value will produce more partial (truncated) terminal segments.\

The recommended value for \*SIZE\* is different for the two datasets:
- For the short dataset, it should be very strict. At least 90% of the reference segment length (e.g. -l 9000 for a 10kb reference segment).
- For the long dataset, it can be more lenient. At least 70% of the reference segment length.

```bash
python3 03_extractSegments.py reads/short.fa.chopped.oriented.fa blastN/reference.vs.short.fa.chop.oriented.sorted.blastN -l *SIZE*
```
```bash
python3 03_extractSegments.py reads/long.fa.chopped.oriented.fa blastN/reference.vs.long.fa.chop.oriented.sorted.blastN -l *SIZE*
```

This script will generate multiple output files.
1. reads/short.fa.chopped.oriented.fa.**report** a quick summary
2. reads/short.fa.chopped.oriented.fa.**segmented.fa** an output fasta file that contains sequences of reads chopped into segments
3. reads/short.fa.chopped.oriented.fa.**undetected.fa** an output fasta file that includes sequences occuring between blast hits in a single read, that are of a minimum size
   - Small gaps between blast hits are common and not neccesarily indicative of a unique sequence (that is not a repeat segment)
   - The parameter --undetectedSegmentMinimumLength can be adjusted to count more or less of these as real features
4. reads/short.fa.chopped.oriented.fa.**allregions.fa** is a combination of file 2 and 3


# 4. Haplotypes

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


# 5. Blocks

Overview:
The most basic ***block*** is a single read listed as the haplotypes of its constituent segments, e.g.

```tsv
read_01  1A 1A 1D 1A 1A
read_03  1B 1D 
read_04  3A
```

In addition, it is very useful to incorporate information about other component sequence (that is not tandem repeat) the read contains.
```tsv
read_01  -- 24 -- --  1A 1A 1D 1A 1A
read_03  -- -- -- --  1B 1D 
read_04  51 -- -- --  3A
```

Now, the first four numbers in each line indicate the percent of the read that contains BLAST hits to four different component sequences.\
In this example, 
1. Left Linker
2. Right Linker
3. Mobile Element
4. Other (not 1., not 2., not 3., and not tandem repeat reference)

Therefore:
- 24 percent of read_01 is Right Linker
- 51 percent of read_04 is Left Linker
- none of the reads contain Mobile Elements
- none of the reads contain Other sequence

The **final blocks file** may look something like this:
```tsv
read_01 l 70243    -- 24 -- --  1A 1A 1D 1A 1A         1B 1D 1F 1B 1B
read_03 l 21015    -- -- -- --  1B 1D                  1G 1f
read_04 l 12972    51 -- -- --  3A                     3C
```
- Column 2 (new) contains the abbreviated id for the dataset
- Column 3 (new) contains the length of the read in base pairs
- Column 4 (new) is empty, allowing for user annotation of block offset in part 6.
- Multiple haplotype lists are shown. For example, for read_01, the first segment is "1A" in haplotype set 1, and "1B" in haplotype set 2.

By gathering all of this information about each read, blocks can be easily aligned together, allowing for assembly of the tandem repeat.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 5a. Assigning haplotypes to each segment

Generate the trace of each segment in the long.fa dataset, using the same methodology as part 4.\
Note that the VAR file here is from variants in the short dataset.

```bash
python3 04_findVariants.py alignments/long.segments.sam alignments/short.segments.sam.VAR.tab reference/reference.fa
```
```bash
bash 04_findTrace.sh alignments/long.segments.sam.variant_table.tab var1
```

Now assign a haplotype to each segment based on the segment's trace, and the haplotypes made in part 4.
```bash
python3 05_classifySegmentHaplotypes.py alignments/long.segments.sam.variant_table.tab.Variants.var1.Traces haplotypes/haplotypes.var1.tsv
```
Do the same for the short reads.
```bash
python3 05_classifySegmentHaplotypes.py alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces haplotypes/haplotypes.var1.tsv
```

The script will output the classification of each segment.
1. If the trace of a segment matches perfectly with a haplotype, it will be assigned that haplotype (soft match)
2. Otherwise, it will be assigned to *the haplotype which trace is the most closest to the segment's* (closest match)

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 5b. Build haplotype lists

```bash
python3 05_combineSegments.py alignments/long.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv
python3 05_combineSegments.py alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv
```

The output file might look something like:
```tsv
68c25dc8-7e91-4130-b22c-07212e7e981e	1O 1o 1o 1l 6d 6D 5a __ 3a 1O 1O 1o 1O 1o 1o 1o 1o 1o 1o 1O 1o 1O 1o 1O 1O 
8dfaf3fb-0886-45a6-b357-121b0fe6ffdc	7j 1p 1P 1P 1P 1P 1i 1p 1P 1p 1P 1P 1P 1P 1P 1P 1p 1P 1P 1p 1P 1p 1P 1P 
8335da6c-1356-4d97-88ea-a4c23f2be420	1o 1O 1l 1L 1O 1o 1l 7E 7M 1k __ 7j 7m 1k __ 7j 7m 7M 7M 7M 7M 7m 7M 
```
Notes:
- The "__" haplotype indicates that there was a significant (100bp or more) undetected region between BLAST hits during segmenting.
- Segments with lowercase letters indicate a "closest match" categorization (see 5a). These could occur because of:
  1. Sequencing errors in the important nucleotides in the trace (most likely)
  2. Real mutations that aren't included in any of the haplotypes (also very plausible)

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 5c. Finding components in reads

The goal here is search for other non-repeat components in the reads.\
This information is indispensable for proper assembly in part 6.\
Any number of components can be chosen to search for in the reads.
It is recommended to prepare at least these sequences:
- left linker: 1000 base pairs on the left end of the tandem repeat cluster
- right linker: 1000 base pairs on the right end of the tandem repeat cluster
In addition, especially in plant genomes, it may be useful to search for the following sequence:
- mobile element coding sequence

Generating references for these sequences is not neccesarily trivial. One may opt to look at contigs that occur at the transition to the tandem repeat cluster.\
If the tandem repeat cluster occurs at the end of a chromosome,  **telomere** can be used for the left or right linker.

```bash
mkdir components
```
Gather the references for these sequences in .fa files (with a single sequence in each):
```bash
cp *path_to_left_linker_fa* components/left_linker.fa
cp *path_to_right_linker_fa* components/right_linker.fa
cp *path_to_mobile_element_fa* components/mobile_element.fa
...
```

For each component:
1. **5c.1** run blastn against the full length reads on each dataset.
```bash
blastn -task blastn -query components/*COMPONENT.fa* -db blastN/DB/short.fa.chop.oriented.BDB -out blastN/*COMPONENT*.vs.short.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```
```bash
blastn -task blastn -query components/*COMPONENT.fa* -db blastN/DB/long.fa.chop.oriented.BDB -out blastN/*COMPONENT*.vs.long.fa.chop.oriented.blastN -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' -evalue 1e-240 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 100 -xdrop_gap 2000 -xdrop_gap_final 4000 -max_target_seqs 64000 -num_threads *THREADS*
```
2. **5c.2** Convert the output to coverage, where each line shows the fraction of that read that is covered by the blast hit.
```bash
bash 05_blast2coverage.sh blastN/*COMPONENT*.vs.short.fa.chop.oriented.blastN reads/short.fa.chop.oriented.fa
```
```bash
bash 05_blast2coverage.sh blastN/*COMPONENT*.vs.long.fa.chop.oriented.blastN reads/long.fa.chop.oriented.fa
```

Finally, convert into coverages for the reference segment:
```bash
bash 05_blast2coverage.sh blastN/reference.vs.short.fa.chop.oriented.blastN reads/short.fa.chop.oriented.fa
```
```bash
bash 05_blast2coverage.sh blastN/reference.vs.long.fa.chop.oriented.blastN reads/long.fa.chop.oriented.fa
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 5d. Combining all information to build blocks


Once there is one .coverage file for each component (and the reference segment), all data can be combined to build **blocks**.

```tsv
bash 05_assembleBlocks.sh NAME \
reference.vs.dataset.coverage \
component1.vs.dataset.coverage \
component2.vs.dataset.coverage \
component3.vs.dataset.coverage \
... \
-- \
dataset.segments.var1.haplotype.list.tab \
dataset.segments.var2.haplotype.list.tab \
...
```
Note:
- There can be any number of components (including 0)
- The -- is used to seperate the components coverage and haplotype list files
- There can be any number of haplotype list files (one for each haplotype set, at least one). The first haplotype list indicates the ***primary haplotype*** set. The second haplotype list is the ***secondary haplotype*** set, etc.

For example, if you have used the components recommended in 5b (left_linker, right_linker, mobile_element), and have two haplotype sets (var1 and var2):
```bash
bash 05_assembleBlocks.sh long \
blastN/reference.vs.long.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/left_linker.vs.long.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/right_linker.vs.long.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/mobile_element.vs.long.fa.chop.oriented.blastN.hits.bed.coverage \
-- \
alignments/long.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv.haplotype.list.tab \
alignments/long.segments.sam.variant_table.tab.Variants.var2.Traces.Classified.tsv.haplotype.list.tab 
```

The output appears in blocks/long.blocks, and should look something like:
```tsv
68c25dc8-7e91-4130-b22c-07212e7e981e	l	244427 		--	--	--	--	1O 1o 1o 1l 6d 6D 5a __ 3a 1O 1O 1o 1O 1o 1o 1o 1o 1o 1o 1O 1o 1O 1o 1O 1O 	1O 1o 1o 1m 5d 5D 5a __ 3a 1o 1O 1o 1O 1o 1o 1o 1o 1o 1o 1O 1o 1O 1o 1O 1O 
8dfaf3fb-0886-45a6-b357-121b0fe6ffdc	l	232951 		--	--	--	--	7j 1p 1P 1P 1P 1P 1i 1p 1P 1p 1P 1P 1P 1P 1P 1P 1p 1P 1P 1p 1P 1p 1P 1P    	1e 1e 1e 1e 1e 1e 1u 1e 1e 1u 1e 1e 1e 1e 1a 1e 1e 1e 1e 1w 1e 1e 1e 1e    
8335da6c-1356-4d97-88ea-a4c23f2be420	l	215403 		--	--	10	--	1o 1O 1l 1L 1O 1o 1l 7E 7M 1k __ 7j 7m 1k __ 7j 7m 7M 7M 7M 7M 7m 7M       	1o 1O 1m 1M 1o 1o 1m 2M 2a 1l __ 3a 2a 1l __ 2k 2a 2a 2a 2a 2a 2a 2a       
6a07f399-43c1-469d-9e3c-a7d2b4401171	l	212265 		--	--	--	--	1J 1P 1P 1P 1P 1P 1P 1P 1P 1P 1P 1P 1P 1P 1p 1P 1P 1p 1p 1P 1P             	1Y 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1c 1Q 1e             
8a45a49e-cdf2-4f76-bc47-e9659b449fd5	l	203958 		--	--	--	--	1P 1P 1P 1p 1P 1P 1P 1P 1P 1p 1P 1p 1p 1P 1A __ 1E 1P 1P 1P 1P             	1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1P __ 1c 1e 1e 1e 1d             
4e2aed6e-afe5-4575-a65d-8df0661c1cc9	l	201124 		--	 0	62	17	7j 7M 7M 7f                                                                	2k 2a 2a 2a                                                                
4ad4da9b-ce9d-438e-a7c7-d72a2a0b5037	l	201084 		--	--	--	--	1P 1P 1P 1P 1P 1P 1P 1p 1p 1p 1p 1P 1P 1P 1p 1P 1p 1P 1P 1A                	1e 1e 1e 1e 1e 1e 1e 1q 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1e 1R                
5b798687-27c9-4356-a3bb-e48a8e7e69f6	l	191318 		--	--	--	--	7M 7m 7E 7m 7m 7M 7h 7e 7e 7g 7m 7e 7M 7m 7M 7m 7M 7M 7E                   	2J 2e 2o 2a 2c 2a 2u 2e 2o 2u 2q 2m 2a 2j 2a 2i 2I 2a 2o                   
ad8533ba-f990-498a-a02e-fefc4b2e5f19	l	189161 		--	--	 2	74	1p 1p 1p 1p                                                                	1e 1e 1e 1e                                                                

```

Notes:
- If the script outputs an empty file, check stdout for errors and double check the names of all of your input files.
- The last column before the haplotype list is the OTHER sequence. This is computed by 1 - (the fraction covered by reference) - (fraction covered by component 1) - (fraction covered by component 2) - ... for each component. Blocks with a value here may contain insertions or link to the rest of the genome. **Alternatively, they may be chimeric reads**.
- For formatting, the default maximum number of segments per block is 25. If the actual number is larger, the variable MOST_SEGMENTS can be modified in 05_assembleBlocks.sh.
- By default, only blocks with three or more segments are retained. This number can be adjusted with the variable MINIMUM_N_SEGMENTS.

Do the same for the short dataset.
```bash
bash 05_assembleBlocks.sh short \
blastN/reference.vs.short.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/left_linker.vs.short.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/right_linker.vs.short.fa.chop.oriented.blastN.hits.bed.coverage \
blastN/mobile_element.vs.short.fa.chop.oriented.blastN.hits.bed.coverage \
-- \
alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv.haplotype.list.tab \
alignments/short.segments.sam.variant_table.tab.Variants.var2.Traces.Classified.tsv.haplotype.list.tab 
```

Finally, combine all datasets into one file.
```bash
cat blocks/long.blocks blocks/short.blocks > blocks/all.blocks
```



# 6. Stacks

Overview: the goal is to manually align, or "stack", as many blocks as possible that contain similar ***motifs***. A motif is a unique sequence of haplotypes, such as **"1i __ 3b"** in the example below.\
Once the entire tandem repeat cluster is stacked, these ***stacked blocks*** will be converted to ***aligned reads*** in part 7. Then, the sequences will be collapsed to generate the consensus draft.\

Condider these blocks that have been stacked (adopted from the assembly of the chromosome 1 NOR region in lettuce, datasets K and O):
```tsv
1   10f081c4-8b1f-4ba5-b0f0-b60a30e5e625	K	89016 	12	--	--	57		-- 3a 3a                                                                      	3b 3d 3d 
2   ef37bf58-b541-4f7a-8b92-bfa989679980	O	44625 	11	--	--	--		-- 3A 3A 3B                                                                   	2f 3D 3D 3f 
3   afeb7684-36d2-4dfd-a27d-db2ea9b730f0	K	58025 	12	--	--	25		-- 3a 3a 3B                                                                   	2f 3d 3d 3f 
4   4590ed0d-a4d2-4511-9af5-0f5ca0d8c0d3	K	68373 	 9	--	--	 2		-- 3a 3a 3b 3b 3b                                                             	2f 3d 3d 3f 3f 3f 
5   cc708727-af1c-4b08-bdae-b7db5f8a6afa	K	49136 	 1	--	--	14		-- 3a 3a 3b                                                                   	2f 3d 3d 3f 
6   5fbb3b74-ba06-452a-95f4-7eca013582c9	K	67035 	 7	--	--	--		-- 3a 3a 3b 3b 3b                                                             	2f 3d 3d 3f 3f 3f 
7   8f83d467-55dd-4ffd-b063-eaee9daaaea5	K	60194 	--	--	--	--		      3a 3b 3b 3b 3a 3B                                                             	3d 3f 3f 3f 3d 3f 
8   5760c0df-37ff-4d9c-854c-c42db042d1d3	K	106266	--	--	--	12 	         3b 3b 3B 3a 3B 3b 3b 3b 3b                                                    	3f 3f 3F 3d 3F 3f 3f 3f 3c 
9   b6cd56a2-9e5a-4bfb-a279-5987db21e9a5	K	70987 	--	--	--	--		            3B 3b 3A 3b 3B 3b 3b                                                          	3f 3f 3d 3f 3f 3f 3f 
10  054e1a09-c1a9-45fe-a6fc-dc8994c9c32a	K	65241 	--	--	--	--		                  3a 3b 3B 3B 3b 3b                                                             	3d 3f 3F 3F 3f 3c 
11  301b2ada-a989-4ab6-9750-7cd67b82253a	O	68905 	--	--	--	17		                     3b 3B 3B 3B 3b                                                                	3f 3F 3F 3F 3C 
12  b0816e37-d95c-4d54-9a95-290803156615	O	106500	--	--	--	13		                  3a 3b 3b 3b 3b 3b -- __ 3b 3b 3b                                              	3d 3f 3f 3c 3f 3c 3a __ 3f 3f 3f 
13  60279406-b117-4390-b415-102e15f05e54	K	105064	--	--	--	13		            3b 3b 3A 3B 3b 3B 3b 3b -- __ 3b                                              	3f 3f 3D 3F 3f 3F 3f 3C 3a __ 3f 
14  ffa8990c-8936-4dbc-8e95-573c3c008f8b	K	161518	--	--	--	58		                                 3b -- __ 3b 3B 3b 3b 3b 3b                                                    	3c 3a __ 3f 3f 3f 3f 3f 3f 
15  20ef918e-ca4f-4084-8928-480bb52f78e6	K	63807 	--	--	--	 8		                                    -- __ 3b 3b 3b 3b 3b                                                          	2f __ 3f 3f 3f 3f 3c 
16  27ef0cf9-2d0c-4078-9785-b1847809e309	K	78385 	--	--	--	18		                        3b 3b 3B 3b -- __ 3b 3B                                                       	3f 3f 3f 3c 3a __ 3f 3F 
17  583acb3f-09a9-4746-a2fe-9636974ef1ff	O	64828 	--	--	--	22		                                 3b -- __ 3b 3B 3B                                                             	3c 3a __ 3f 3F 3F 
18  f6de92c8-744f-4c83-8595-2499d01634b6	K	69377 	--	--	--	20		                        3b 3B 3B 3b -- __ 3b                                                          	3f 3F 3f 3C 3a __ 3f 
19  43e8b987-6b98-4657-ad70-8cd5ffd2f5d1	O	172846	--	--	--	12		                              3B 3b 1i __ 3b 3B 3b 3B 3B 3b __ 2C 4a 4a 4A 4C 4c 4c 4c                      	3f 3C 3a __ 3f 3F 3f 3F 3F 3f __ 3A 4b 4A 4B 4K 4k 4f 4h 
20  995458f7-f1f1-4486-b4c8-686d5a9a1dbd	O	168857	--	--	--	13		                           3B 3B 3b 1i __ 3b 3B 3b 3b 3B 3b __ 2C 4A 4a 4C 4C 4c 4c                      	3F 3F 3C 3a __ 3f 3F 3f 3f 3F 3f __ 3A 4B 4a 4K 4K 4k 4i 
21  fc9b176f-0a65-41d4-b05f-17f8c50afca0	O	137969	--	--	--	16		                              3B 3b 1i __ 3b 3B 3B 3B 3B 3b __ 2C 4a 4a 4C                                  	3F 3c 3a __ 3f 3F 3F 3F 3F 3f __ 3A 4b 4a 4K 
22  45deef91-e9ea-4810-98f0-52573132833c	O	127339	--	--	--	17		                        3b 3b 3B 3b 1i __ 3b 3b 3B 3b 3b 3b                                           	3f 3f 3F 3c 3a __ 3f 3f 3F 3f 3f 3f 
23  c7798e2f-0d6b-49e8-a09c-6401fdde1e4a	O	111917	--	--	--	12		            3b 3B 3A 3B 3b 3b 3B 3b 1i __ 3b 3b                                           	3f 3F 3d 3F 3f 3f 3f 3C 3a __ 3f 3f 
```
- The first 6 blocks contain left linker sequence, indicating this is the beginning of the tandem repeat cluster
- The 7th block, contains the motif "3a 3b 3b 3b 3a" which is used to expand the Traces
- Block 12 connects the "3a" motif and reveals a new motif, "-- __ 3b", which can also appear as "1i __ 3b"
- Block 19 connects the "-- __ 3b" motif to a new motif, "3b __ 2C", both which are retrotransposon insertions
- Blocks 12 through 23 all contain "other" sequence, indicating the sequence in retrotransposon insertion occuring at the undetected segment "__"
- The secondary haplotypes also stack perfectly. 

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6a. Preprocess blocks file


By default, the blocks are sorted from longest to shortest. It may be useful to float blocks containing left_linker and/or right_linker (fields 4 and 5 in the example above) to the top (or any other components):
```bash
sort -k4,4nr -k5,5nr -k3,3nr blocks/all.blocks >  blocks/all.sorted.blocks
```
In addition, it may be quite useful to group the blocks by which haplotype groups they contain.
```bash
python3 06_groupBlocks.py blocks/all.sorted.blocks
```
<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6b. Stacking blocks

```bash
cp blocks/all.sorted.grouped.blocks > blocks/all.stacks
echo "" > blocks/all.unused.stacks
```

This is the bulk of the work in this workflow.\
**Manually add space characters directly before the primary haplotype list in *all.sorted.grouped.stacks* to align their motifs together.**
- Stacking should start with reads containing the left linker (the start of the tandem repeat cluster), then expand outwards until reaching the right linker (the end of the tandem repeat cluster). 
- **Each segment should have at least 10 depth covering it**, with virtually identical haplotypes in all haplotype sets.
  - In the example above, the first few segments (-- 3A 3A) have only 6, 6, and 7 depth (respectively). More blocks need to be stacked.
  - Especially with high sequencing coverage, not every block needs to be stacked. Unstacked blocks should be placed into *blocks/all.sorted.grouped.unused.blocks*.
    - *Any blocks that are moved here should be verified that the motifs they contain already exist in the assembled block stacks*.
- Visual Studio code is a great text editor for this task.
  - Tab Size 3
  - RaindbowCSV Extension
    - Modify the field colors to make the field with primary haplotypes stand out
- The segments at the ends of blocks should always be used with caution. They may be truncated and their haplotypes may be innacurate and misleading.
- Try to find motifs and string them into longer regions
- Once you have many large regions, they can be connected together
- In many cases, the quantity and exact nature of the connection between these large regions may not be known. In this case, there are a few options:
   1. it may be practical to fill these gaps with consensus segments (see below)
   2. these large regions can be kept seperate and assembled as seperate contigs

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6c. Consensus segments

Tandem repeat clusters may have long regions of repeats that are practically identical to the point of indistinguishability in haplotypes.\
> *These sequences may contain active genes and/or are under selection*

These sequences can be assembled by manually stacking blocks, but this can be tedious.\
Instead, a consensus sequence can be used to automatically fill in multiple segments in a row with a copy and paste sequence.

Consider this example (adopted from the assembly of the chromosome 1 NOR region in lettuce)

```tsv
1    X-REPEAT_LEN----------------------------	-	9847
2    X-SEGMENTLEN----------------------------	-	7000
3    X-COMPONENTS----------------------------	-	telomere,copia_retrotransposon,right_linker
4    X-DEFINE--------------------------------	2K	PATH_TO_FOLDER/reference/2K.sam
5    
6    #### Contig 1
7    # Adopted from assembly of chromosome 1 NOR in lettuce
8
9    a11c9772-e066-4d14-b793-81f457c8d34f    	K	64056 	--	--	--	11		         2K 2e 2k 2k 1i __ 2i                                                          	2Z 2e 2z 2z 1t __ 2f 
10   0f479704-621d-445e-81eb-1c04ecfe4936    	K	64072 	--	--	--	11		         2k 2k 2i 2k 1i __ 2i                                                          	2z 2z 2q 2z 1t __ 2f 
11   ada4bfc0-a05a-45cb-ba53-171eb085d974    	K	64497 	--	--	--	11		         2k 2k 2K 2e 1i __ 2i 2f                                                       	2z 2z 2Z 2e 1t __ 2f 2y 
12   cc8fcd4b-6580-4f0c-81ca-9a2099de66db    	K	67472 	--	--	--	11		               2K 2k 1i __ 2i 2K __ 2k 2k                                                    	2Z 2k 1t __ 2f 2y __ 2x 2x 
13   e8a1af28-d84b-4011-bf8b-e0f161fabb9d    	K	55934 	--	--	--	13		            2K 2K 2D 1i __ 2i 2K                                                          	2Z 2Z 2D 1t __ 2f 2Y 
14   X-INSERT--------------------------------	2K	00002 	--	--	--	--		                              05 05
15   X-INSERT--------------------------------	2K	00020 	--	--	--	--		                                    10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 
16   X-INSERT--------------------------------	2K	00003 	--	--	--	--		                                                                                                05 05 05
17   4e2aed6e-afe5-4575-a65d-8df0661c1cc9    	O	201124	--	 7	 0	72		                                                                                                   2i 2K 2K 2h                                                                   	2o 2Z 2Z 2j 
18   188162e8-a68c-401d-af9a-87c2c5ba33c3    	O	133120	--	--	 0	 4		                                                                           2K 2k 2k 2K 2K 2k_2k 2k 2k 2k 2k 2K                                           	2Z 2z 2z 2Z 2z 2z_2z 2z 2z 2p 2z 2Z 
19   874896aa-9147-414c-8517-fa4d34023f2b    	O	61776 	--	--	 2	33		                                                                                                   -- 2K 2K 2K                                                                   	2f 2Z 2Z 2Z 
20
21   #### Contig 2
22   # Example of creating another (very small) contig from another collection of stacked blocks
23
24   50f67900-7159-405e-beea-3cd4943d3a3f    	K	155823	--	--	--	--		1n 1n 1N 1N 1N 1N 1N 1N 1N 1N 1N 1I __ 2C 1n 1F 1N                            	1z 1z 1Z 1z 1Z 1Z 1e 1Z 1Z 1A 1Z 1T __ 3A 1z 1O 1Z 
25   77974855-4275-46c0-b4e7-ba0d39ae3359    	O	140381	--	--	--	--		            1n 1n 1N 1N 1N 1N 1N 1I __ 2C 1N 1F 1e 1f 1N 1F 1n                            	1N 1z 1Z 1Z 1Z 1Z 1Z 1T __ 3A 1Z 1O 1a 1o 1Z 1O 1a 
26   90ea782e-718c-49ad-ae96-eb8319f63d30    	K	127630	--	--	--	--		                        1n 1N 1N 1i __ -- 1n 1f 1N 1f 1n 1F 1n 1n 1n                                  	1z 1Z 1Z 1t __ 3a 1z 1o 1Z 1o 1z 1O 1z 1z 1a 
27   5f9c5373-ece6-4fcf-b9db-f034b923e864    	K	127521	--	--	--	--		                        1n 1N 1N 1I __ 2C 1C 1F 1N 1f 1n 1f 1N 1N 1n                                  	1z 1Z 1z 1T __ 3A 1M 1o 1Z 1o 1e 1o 1Z 1z 1a 
28   fcd3506c-1756-420a-9bf3-fd4ade2a17b8    	O	120105	--	--	--	--		                              1N 1I __ 2C 1N 1F 1N 1N 1N 1f 1N 1f 1n 1n                                     	1Z 1T __ 3A 1z 1O 1Z 1Z 1Z 1o 1Z 1o 1z 1z 
```

- It was determined that there was a repition of the "2K" haplotype an unknown number of times between the motif shown on lines 6-9 and the ending of the tandem repeat cluster (the blocks on lines 14 - 16).

- Lines 1, 2 and 3 are required headers that define the length of the repeat, minimum size of terminal segments, and the components fields (see 6d).

- On line 1, the consensus sequence is defined. This special "block" contains the following fields:
    1. "X-DEFINE----------------------------" indicates to the parser that this is defining a consensus sequence
    2. The name of the consensus sequence. In this case, it represents the common 2K haplotype sequence.
    3. An absolute file path to the consensus sequence aligned to the reference segment:
      ```bash
      cp *PATH_TO_CONSENSUS* /referece/2K.consensus.fa
      minimap2 --eqx -a reference/reference.fa reference/2K.fa > reference/2K.sam
      ```

- Line 3 indicates the start of a contig (see 6d)

- On line 11, the consensus sequence is inserted. This special "block" contains the following fields:
    1. "X-INSERT----------------------------" indicates to the parser that this is inserting a repeated consensus sequence
    2. The name of the consensus sequence. It should match one of the previously defined names.
    3. The number of segments to insert in a row. For line 11, it is two segments in a row.
    4-7. Empty
    8. The **location** in which to insert the consensus segments is determined by where the number begins. The **depth** of each insertion (more information in part 7) if determined by the value of the first number. It is not neccesary to copy and paste this number , but doing so helps visually indicate where the insertion ends.

- Lines 12 and 13 once again define consensus insertions. There is a total of 25 segments that have consensus insertions.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6d. Defining contigs and adding headers

Add the following three required headers (tab seperated) to the beginning of the stacked blocks file:
```tsv
X-REPEAT_LEN------------------------	-	9847
X-SEGMENTLEN------------------------	-	7000
X-COMPONENTS------------------------	-	telomere,copia_retrotransposon,right_linker
```
1. X-REPEAT_LEN - The length of the reference segment (here 9847 base pairs)
2. X-SEGMENTLEN - The precise value used from step3b to determine the minimum length of a segment for the *long* dataset. This value is important to help determine if the first segment in each block is a partial segment, or a full segment (and the partial segment was removed at step 3b).
3. X-COMPONENTS - A comma seperate list of components used in 5d.

If the entire tandem repeat cluster couldn't be assembled in a single large region, it may be split into multiple contigs. Each contig should be defined by a header starting with "#### ", like on lines 3 and 18 above. All blocks preceding that header (until the next header) will belong to that contig.\
```tsv
#### Contig_Name
```

Comments may be addded (starting with "# ") to describe different locations.
```tsv
# Comment
```

Once the blocks/all.stacks file is finalized, copy it over. V1 can be replaced with V2, etc. for future versions.
```bash
cp blocks/all.stacks >  blocks/all.V1.stacks
```

# 7. Draft

Overview: generate a draft of the tandem repeat cluster DNA sequence by collapsing SAM alignments that are mapped by parsing where the blocks have been manually aligned in part 6.

```bash
mkdir drafts
```
<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 7a. Generate SAM alignments onto tandem reference sequence.

First, tandem duplicate the reference sequence.\
Apprxoimate the total segments that the longest single contig in *blocks/all.V1.stacks* has, then add 10 more.\
Replace \*TOTAL_COUNT\* with this value.
```bash
awk -v VAR1=*TOTAL_COUNT* 'NR == 1 {print $0} NR == 2 {for(i = 0; i < VAR1; i++) {printf($0)};}' reference/reference.fa > reference/reference.tandem.fa
```
Duplicate the annotations file. Replace \*REFERENCE_LEN\* with the length of your reference.
```bash
awk -v VAR1=*TOTAL_COUNT* -v VAR2=*REFERENCE_LEN* 'NR == 1 {print $0; next} {for(i = 0; i < *TOTAL_COUNT*; i++) {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%d\t\n",$1,$2,$3,$4 + VAR2*i,$5 + VAR2*i,$6,$7,$8,$9,i+1)};}' reference/annotations.gff > reference/annotations.tandem.gff
```

Combine all full length read datasets.
```bash
cat reads/long.fa.chop.oriented.fa reads/short.fa.chop.oriented.fa > reads/all.fa
```
Align the full length reads from all datasets to the duplicated reference sequence.
```bash
minimap2  -t *THREADS* --eqx -a -N 1 -E 2,0 --end-bonus=8 reference/reference.tandem.fa reads/all.fa > drafts/all.tandem.sam
```
Retain only primary alignments.
```bash
awk 'and($2, 2304) == 0 {print}' drafts/all.tandem.sam > drafts/all.tandem.primary.sam
```

These reads are being mapped to a tandem repeat, so their mapping position is arbitrary and chosen randomly.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 7b. Assign mapping positions from stacked blocks file.

```bash
python3 07_stackedBlocksToSAM.py blocks/all.V1.stacks drafts/all.tandem.primary.sam
```

For each contig header *contig_name*, one output file will be produced drafts/*contig_name*.sam.\
This file contains one read for each block in the input file, with adjusted mapping positions.\

Visualize each SAM file. The visualization could look something like this:

![7b_output_example1](/assets/07b_example1.png)
> *Start of Lettuce Chromosome 8 NOR region*. Visualization with CLC Genomics Workbench. Polymorphisms with respect to tandem duplicated reference indicated by colored bars. Note soft clipping of left-linker sequence on left border. Blocks are stacked with up to 58 depth.
![7b_output_example2](/assets/07b_example2.png)
> *Middle of Lettuce Chromosome 8 NOR region*. Visualization with CLC Genomics Workbench. Large insertion of mobile element.


