# 5. Constructing Blocks

Overview:
The most basic ***block*** is a single read listed as the haplotypes of its constituent segments, e.g.

```tsv
read_01  1A 1A 1D 1A 1A
read_03  1B 1D 
read_04  3A
```

In addition, it is very useful to incorporate information about other component sequences (that are not tandem repeat) the read contains.
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
python3 scripts/04_findVariants.py alignments/long.segments.sam alignments/short.segments.sam.VAR.tab reference/reference.fa
```
```bash
bash scripts/04_findTrace.sh alignments/long.segments.sam.variant_table.tab var1
```

Now assign a haplotype to each segment based on the segment's trace, and the haplotypes made in part 4.
```bash
python3 scripts/05_classifySegmentHaplotypes.py alignments/long.segments.sam.variant_table.tab.Variants.var1.Traces haplotypes/haplotypes.var1.tsv
```
Do the same for the short reads.
```bash
python3 scripts/05_classifySegmentHaplotypes.py alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces haplotypes/haplotypes.var1.tsv
```

The script will output the classification of each segment.
1. If the trace of a segment matches perfectly with a haplotype, it will be assigned that haplotype (soft match)
2. Otherwise, it will be assigned to *the haplotype which trace is the most closest to the segment's* (closest match)

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 5b. Build haplotype lists

```bash
python3 scripts/05_combineSegments.py alignments/long.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv
python3 scripts/05_combineSegments.py alignments/short.segments.sam.variant_table.tab.Variants.var1.Traces.Classified.tsv
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

Generating references for these sequences is not necessarily trivial. One may opt to look at contigs that occur at the transition to the tandem repeat cluster.\
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
bash scripts/05_blast2coverage.sh blastN/*COMPONENT*.vs.short.fa.chop.oriented.blastN reads/short.fa.chop.oriented.fa
```
```bash
bash scripts/05_blast2coverage.sh blastN/*COMPONENT*.vs.long.fa.chop.oriented.blastN reads/long.fa.chop.oriented.fa
```

Finally, convert into coverages for the reference segment:
```bash
bash scripts/05_blast2coverage.sh blastN/reference.vs.short.fa.chop.oriented.blastN reads/short.fa.chop.oriented.fa
```
```bash
bash scripts/05_blast2coverage.sh blastN/reference.vs.long.fa.chop.oriented.blastN reads/long.fa.chop.oriented.fa
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
bash scripts/05_assembleBlocks.sh long \
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
bash scripts/05_assembleBlocks.sh short \
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

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-06-STACKS](WORKFLOW-06-STACKS.md)