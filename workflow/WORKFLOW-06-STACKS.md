# 6. Stacking Blocks

Overview: the goal is to manually align, or "stack", as many blocks as possible that contain similar ***motifs***. A motif is a unique sequence of haplotypes, such as **"1i __ 3b"** in the example below.\
Once the entire tandem repeat cluster is stacked, these ***stacked blocks*** will be converted to ***aligned reads*** in part 7. Then, the sequences will be collapsed to generate the consensus draft.\

Consider these blocks that have been stacked (adopted from the assembly of the chromosome 1 NOR region in lettuce, datasets K and O):
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
- Blocks 12 through 23 all contain "other" sequence, indicating the sequence in retrotransposon insertion occurring at the undetected segment "__"
- The secondary haplotypes also stack perfectly. 

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6a. Preprocess blocks file


By default, the blocks are sorted from longest to shortest. It may be useful to float blocks containing left_linker and/or right_linker (fields 4 and 5 in the example above) to the top (or any other components):
```bash
sort -k4,4nr -k5,5nr -k3,3nr blocks/all.blocks >  blocks/all.sorted.blocks
```
In addition, it may be quite useful to group the blocks by which haplotype groups they contain.
```bash
python3 scripts/06_groupBlocks.py blocks/all.sorted.blocks
```
<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6b. Stacking blocks

```bash
cp blocks/all.sorted.grouped.blocks blocks/all.stacks
echo "" > blocks/all.unused.stacks
```

This is the bulk of the work in this workflow.\
**Manually add space characters directly before the primary haplotype list in *all.stacks* to align their motifs together.**
- Stacking should start with reads containing the left linker (the start of the tandem repeat cluster), then expand outwards until reaching the right linker (the end of the tandem repeat cluster). 
- **Each segment should have at least 10 depth covering it**, with virtually identical haplotypes in all haplotype sets.
  - In the example above, the first few segments (-- 3A 3A) have only 6, 6, and 7 depth (respectively). More blocks need to be stacked.
  - Especially with high sequencing coverage, not every block needs to be stacked. Unstacked blocks should be placed into *blocks/all.sorted.grouped.unused.blocks*.
    - *Any blocks that are moved here should be verified that the motifs they contain already exist in the assembled block stacks*.
- Each block should appear exactly once between the files *all.stacks* and *all.unused.stacks*.
- Visual Studio Code text editor is highly recommended for this step.
  - Set Tab Size to 3
  - Install RainbowCSV Extension
    - Modify the field colors to make the field with primary haplotypes stand out
- The segments at the ends of blocks should always be used with caution. They may be truncated and their haplotypes may be inacurate and misleading.
- Try to find motifs and string them into longer regions
- Once you have many large regions, they can be connected together
- In many cases, the quantity and exact nature of the connection between these large regions may not be known. In this case, there are a few options:
   1. It may be practical to fill these gaps with consensus segments (see below).
   2. These large regions can be kept separate and assembled as separate contigs.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6c. Consensus segments

Tandem repeat clusters may have long regions of repeats that are practically identical to the point of indistinguishability in haplotypes.\
> *These sequences may contain active genes and/or are under selection*

These sequences can be assembled by manually stacking blocks, but this can be tedious.\
Instead, a consensus sequence can be used to automatically fill in multiple segments in a row with a copy and paste sequence.

Consider this example (adopted from the assembly of the chromosome 1 NOR region in lettuce)

```tsv
     1	X-REPEAT_LEN----------------------------	-	9847
     2	X-SEGMENTLEN----------------------------	-	7000
     3	X-COMPONENTS----------------------------	-	telomere,copia_retrotransposon,right_linker
     4	X-DEFINE--------------------------------	7M	/share/rwmwork/peter/Ribo_Variants/test/reference/7M.sam
     5	
     6	#### Contig_1
     7	# Adopted from assembly of chromosome 8 NOR in lettuce
     8	
     9	a11c9772-e066-4d14-b793-81f457c8d34f    	l	64056  	--	--	11	--		   7M 7f 7m 7m 1k __ 7j                                                       	2a 2a 2a 2j 1l __ 2x                                                       
    10	0f479704-621d-445e-81eb-1c04ecfe4936    	l	64072  	--	--	11	--		   7m 7m 7k 7m 1k __ 7j                                                       	2a 2h 2k 2a 1l __ 3a                                                       
    11	ada4bfc0-a05a-45cb-ba53-171eb085d974    	l	64497  	--	--	11	--		   7A 7m 7M 7m 1k __ 7J 7m                                                    	2a 2a 2a 2a 1l __ 3a 2x                                                    
    12	cc8fcd4b-6580-4f0c-81ca-9a2099de66db    	l	67472  	--	--	11	--		         7M 7m 1k __ 7J 7M __ 7m 7l                                                 	2a 2a 1l __ 3a 2x __ 2a 2z                                                 
    13	e8a1af28-d84b-4011-bf8b-e0f161fabb9d    	l	55934  	--	--	13	--		      7M 7M 7m 1k __ 7j 7M                                                       	2a 2a 2a 1l __ 3a 2x                                                       
    14	X-INSERT--------------------------------	7M	00002 	--	--	--	--		                        05 05
    15	X-INSERT--------------------------------	7M	00020 	--	--	--	--		                              10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 
    16	X-INSERT--------------------------------	7M	00003 	--	--	--	--		                                                                                          05 05 05
    17	4e2aed6e-afe5-4575-a65d-8df0661c1cc9    	l	201124 	--	 0	62	17		                                                                                             7j 7M 7M 7f                                                                	2k 2a 2a 2a                                                                
    18	188162e8-a68c-401d-af9a-87c2c5ba33c3    	l	133120 	--	 0	 0	 3		                                                                     7M 7m 7m 7M 7M 7m 7m 7m 7m 7m 7m 7M                                        	2H 2a 2a 2a 2a 2a 2a 2a 2a 2a 2a 2a                                        
    19	874896aa-9147-414c-8517-fa4d34023f2b    	l	61776  	--	 2	25	 8		                                                                                                7m 7M 7M                                                                   	2a 2a 2a                                                                   
    20	
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

- It was determined that there was a repetition of the "2K" haplotype an unknown number of times between the motif shown on lines 6-9 and the ending of the tandem repeat cluster (the blocks on lines 14 - 16).

- Lines 1, 2 and 3 are required headers that define the length of the repeat, minimum size of terminal segments, and the components fields (see 6d).

- On line 4, the consensus sequence is defined. This special "block" contains the following fields:
    1. "X-DEFINE----------------------------" indicates to the parser that this is defining a consensus sequence
    2. The name of the consensus sequence. In this case, it represents the common 2K haplotype sequence.
    3. An absolute file path to the consensus sequence aligned to the reference segment:
      ```bash
      cp *PATH_TO_CONSENSUS* /referece/2K.consensus.fa
      minimap2 --eqx -a reference/reference.fa reference/2K.fa > reference/2K.sam
      ```

- Line 3 indicates the start of a contig (see 6d)

- On line 14, the consensus sequence is inserted. This special "block" contains the following fields:
    1. "X-INSERT----------------------------" indicates to the parser that this is inserting a repeated consensus sequence
    2. The name of the consensus sequence. It should match one of the previously defined names.
    3. The number of segments to insert in a row. For line 14, it is two segments in a row.
    4-7. Empty
    8. The **location** in which to insert the consensus segments is determined by where the number begins. The **depth** of each insertion (more information in part 7) if determined by the value of the first number. It is not necessary to copy and paste this number , but doing so helps visually indicate where the insertion ends.

- Lines 15 and 16 once again define consensus insertions. There are a total of 25 segments that have consensus insertions.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 6d. Defining contigs and adding headers

Add the following three required headers (tab separated) to the beginning of the stacked blocks file:
```tsv
X-REPEAT_LEN------------------------	-	9847
X-SEGMENTLEN------------------------	-	7000
X-COMPONENTS------------------------	-	telomere,copia_retrotransposon,right_linker
```
1. X-REPEAT_LEN - The length of the reference segment (here 9847 base pairs)
2. X-SEGMENTLEN - The precise value used from step 3b to determine the minimum length of a segment for the *long* dataset. This value is important to help determine if the first segment in each block is a partial segment, or a full segment (and the partial segment was removed at step 3b).
3. X-COMPONENTS - A comma separated list of components used in 5d.

If the entire tandem repeat cluster couldn't be assembled in a single large region, it may be split into multiple contigs. Each contig should be defined by a header starting with "#### ", like on lines 6 and 21 above. All blocks preceding that header (until the next header) will belong to that contig.\
```tsv
#### Contig_Name
```

Comments may be added (starting with "# ") to describe different regions (they will have no effect on output).
```tsv
# Comment
```

Once the blocks/all.stacks file is finalized, copy it over. V1 can be replaced with V2, etc. for future versions.
```bash
cp blocks/all.stacks blocks/all.V1.stacks
```

Next step: ![WORKFLOW-07-DRAFT](WORKFLOW-07-DRAFT.md)