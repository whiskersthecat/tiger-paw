# 2. Process Reads

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
python3 scripts/02_chopDuplexAndOrientReads.py reads/short.fa blastN/reference.vs.short.sorted.blastN --include_original
python3 scripts/02_chopDuplexAndOrientReads.py reads/long.fa blastN/reference.vs.long.sorted.blastN --include_original
```

This script will generate multiple output files.
1. reads/short.fa.**report** a count of how many reads were in each of 5 categories
2. reads/short.fa.**readids** the categorization of each read
3. reads/short.fa.**chop.oriented.fa** an output fasta file that includes modified reads:
   - With --include_original, any reads marked as duplexes will have two corresponding reads in the output file:
      - read_name: the modified read.
      - read_name_dp (duplex perfect) or read_name_dbf (duplex bias forward) or read_name_dbr (duplex bias reverse). This read has the sequence of the original read, and indicates its duplex categorization.
   - Otherwise, only modified reads are output

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-03-SEGMENTS](WORKFLOW-03-SEGMENTS.md)