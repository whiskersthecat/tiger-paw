# 7. Creating Draft

Overview: generate a draft of the tandem repeat cluster DNA sequence by collapsing SAM alignments that are mapped by parsing where the blocks have been manually aligned in part 6.

```bash
mkdir drafts
```
<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 7a. Generate SAM alignments onto tandem reference sequence.

First, tandem duplicate the reference sequence.\
Approximate the total segments that the longest single contig in *blocks/all.V1.stacks* has, then add 10 more (this number does not need to be exact).\
Replace \*TOTAL_COUNT\* with this value.
```bash
awk -v VAR1=*TOTAL_COUNT* 'NR == 1 {print $0} NR == 2 {for(i = 0; i < VAR1; i++) {printf($0)};}' reference/reference.fa > reference/reference.tandem.fa
```
Duplicate the annotations file. Replace \*REFERENCE_LEN\* with the length of your reference.
```bash
awk -v VAR1=*TOTAL_COUNT* -v VAR2=*REFERENCE_LEN* 'NR == 1 {print $0; next} {for(i = 0; i < VAR1; i++) {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;REPEAT=%d\t\n",$1,$2,$3,$4 + VAR2*i,$5 + VAR2*i,$6,$7,$8,$9,i+1)};}' reference/features.gff > reference/features.tandem.gff
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

![7b_output_example1](/assets/07b_example2.png)
> *Start of Lettuce Chromosome 8 NOR region*. Visualization with CLC Genomics Workbench. Polymorphisms with respect to tandem duplicated reference indicated by colored bars. Note soft clipping of left-linker sequence on left border. Blocks are stacked with up to 58 depth.

![7b_output_example2](/assets/07b_example1.png)
> *Middle of Lettuce Chromosome 8 NOR region*. Visualization with CLC Genomics Workbench. Large insertion of mobile element.



<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 7c. Tweak mapping positions after large deletion / insertion motifs.

Large insertions and deletions may offset mapping positions in downstream blocks.\
Consider these stacked blocks:

```tsv
    21	#### Contig_2
    22	# Adopted from assembly of chromosome 1 NOR in lettuce
    23	
    24	acdd29f0-6088-4f15-bde1-3b5feeee8fa0    	l	94020  	--	--	15	--		   4C 4C 4C 4c 4c 1k __ 4c 4c                                                 	4A 4a 4A 4a 4a 3a __ 4a 4a                                                 
    25	45deef91-e9ea-4810-98f0-52573132833c    	l	127339 	--	--	17	--		      4c 4c 4C 4c 1k __ 4c 4c 4C 4c 4c 4b                                        	4a 4a 4A 4a 3a __ 4a 4a 4A 4a 4a 4a                                        
    26	43e8b987-6b98-4657-ad70-8cd5ffd2f5d1    	l	172846 	--	--	12	--		            4C 4c 1k __ 4c 4C 4c 4C 4C 4B __ 3A 5a 5a 5A 6D 6d 6d 6d                   	4a 4a 3a __ 4a 4A 4a 4a 4a 4a __ 3A 5a 5a 5A 5D 5d 5d 5d                   
    27	6a79d732-37e9-4441-ae71-f38714ba7bf7    	l	119802 	--	--	18	--		                  1k __ 4c 4c 4c 4c 4c 4b __ 3A 5a 5a 6d                                     	3a __ 4a 4a 4a 4a 4a 4a __ 3A 5a 5a 5d                                     
    28	995458f7-f1f1-4486-b4c8-686d5a9a1dbd    	l	168857 	--	--	13	--		         4C 4C 4c 1k __ 4c 4C 4c 4c 4C 4b __ 3A 5A 5a 6D 6D 6d                      	4A 4A 4a 3a __ 4a 4A 4a 4a 4A 4a __ 3A 5A 5a 5D 5D 5d                      
    29	fc9b176f-0a65-41d4-b05f-17f8c50afca0    	l	137969 	--	--	16	--		            4C 4c 1k __ 4c 4C 4C 4C 4C 4B __ 3A 5a 5a 6D                               	4A 4a 3a __ 4a 4A 4A 4A 4A 4a __ 3A 5a 5a 5D                               
    30	89e3f606-9c33-45d1-bf13-501e6ccb3a83    	l	99200  	--	--	 7	--		                                 4c 4c 4b __ 3A 5a 5a 6d 6d 6d 6d                                           	4a 4a 4a __ 3A 5a 5a 5d 5d 5d 5d                                           
    31	9bfd6db4-35a4-41c1-b963-b224652ce643    	l	95009  	--	--	 8	--		                              4c 4c 4c 4b __ 3A 5a 5a 6D 6D                                              	4a 4a 4a 4a __ 3A 5a 5a 5d 5d                                              
```

The blocks appear to be correctly aligned. As it turns out, the mapping positions blocks on line 30 and 31 are actually off by two segments:

![7b_output_example3](/assets/07b_example3.png)
> For the "4b __ 3A" motif (a large mobile element insertion), most of the blocks are stacking correctly in the Contig_2.sam file. However, two of these blocks (lines 30 and 31) are mapping two segments to the right.

This occurs because of original the motif "1k __ 4c":
- While this motif consumes 3 segments as a block, it only represents 1 actual segment in the alignment of the read.
- The "__" segment is the undetected segment, which relates to a 14,000 base pair insertion.
- The "1k" and "4c" segments come from the truncated segments to the left and right of the insertion that were found when the read was chopped into segments.
- In order to correctly map any read that that is positioned past this segment, its **shift** field for that block must be set to -2. This value tells the parser to shift those blocks left to segments when determining the mapping position.
  - This field allows the blocks to appear stacked (visually) in the stacked blocks file, while mapping correctly when converted to SAM.

Consider this adjustment (note that -2 on lines 30 and 31):
```tsv
    21	#### Contig_2
    22	# Adopted from assembly of chromosome 1 NOR in lettuce
    23	
    24	acdd29f0-6088-4f15-bde1-3b5feeee8fa0    	l	94020  	--	--	15	--		   4C 4C 4C 4c 4c 1k __ 4c 4c                                                 	4A 4a 4A 4a 4a 3a __ 4a 4a                                                 
    25	45deef91-e9ea-4810-98f0-52573132833c    	l	127339 	--	--	17	--		      4c 4c 4C 4c 1k __ 4c 4c 4C 4c 4c 4b                                        	4a 4a 4A 4a 3a __ 4a 4a 4A 4a 4a 4a                                        
    26	43e8b987-6b98-4657-ad70-8cd5ffd2f5d1    	l	172846 	--	--	12	--		            4C 4c 1k __ 4c 4C 4c 4C 4C 4B __ 3A 5a 5a 5A 6D 6d 6d 6d                   	4a 4a 3a __ 4a 4A 4a 4a 4a 4a __ 3A 5a 5a 5A 5D 5d 5d 5d                   
    27	6a79d732-37e9-4441-ae71-f38714ba7bf7    	l	119802 	--	--	18	--		                  1k __ 4c 4c 4c 4c 4c 4b __ 3A 5a 5a 6d                                     	3a __ 4a 4a 4a 4a 4a 4a __ 3A 5a 5a 5d                                     
    28	995458f7-f1f1-4486-b4c8-686d5a9a1dbd    	l	168857 	--	--	13	--		         4C 4C 4c 1k __ 4c 4C 4c 4c 4C 4b __ 3A 5A 5a 6D 6D 6d                      	4A 4A 4a 3a __ 4a 4A 4a 4a 4A 4a __ 3A 5A 5a 5D 5D 5d                      
    29	fc9b176f-0a65-41d4-b05f-17f8c50afca0    	l	137969 	--	--	16	--		            4C 4c 1k __ 4c 4C 4C 4C 4C 4B __ 3A 5a 5a 6D                               	4A 4a 3a __ 4a 4A 4A 4A 4A 4a __ 3A 5a 5a 5D                               
    30	89e3f606-9c33-45d1-bf13-501e6ccb3a83    	l	99200  	--	--	 7	--	-2	                                 4c 4c 4b __ 3A 5a 5a 6d 6d 6d 6d                                           	4a 4a 4a __ 3A 5a 5a 5d 5d 5d 5d                                           
    31	9bfd6db4-35a4-41c1-b963-b224652ce643    	l	95009  	--	--	 8	--	-2	                              4c 4c 4c 4b __ 3A 5a 5a 6D 6D                                              	4a 4a 4a 4a __ 3A 5a 5a 5d 5d                                              
```
Now, the output SAM file looks like this:
![7b_output_example4](/assets/07b_example4.png)
> Now segments 30 and 31 are correctly stacking! 
- Any other blocks that are aligned any position past the "1k _ 4c" will need to have a shift of -2.
- The motif "4b __ 3A" has a similar case. In fact, any blocks that are aligned past this position will need to have a shift of -4.

**Go back to step 6 and adjust the shifts of blocks until they are all mapping correctly**. Then procede to 7d.

Other tips:
- Some large deletions may also consume multiple segments as blocks.
  - Some may consume 2, while others may consume 3.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 7d. Collapsing alignments to generate consensus.

Now that reads have been mapped to their correct locations on the tandem repeat, a consensus can be formed.\
First, call variants on the SAM file:

> TO BE COMPLETED

