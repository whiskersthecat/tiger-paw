# 3. Create Segments

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

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-04-HAPLOTYPES](WORKFLOW-04-HAPLOTYPES.md)