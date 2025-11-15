# 9. Incorporation

Overview: incorporate the assembled tandem repeat cluster contigs(s) into the chromosomal assembly.\
Original chromosome-scale assembly:
```tsv
chrleft - GAP - chrright
~ OR ~
chrleft - COLLAPSEDTRC - chrright
```

New chromosome-scale assembly (fixed):
```tsv
chrleft - TRC - chrright
```

```bash
mkdir genome
```
<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 9a. Prepare sequences.

For each tandem repeat cluster contig generated in step 8, prepare the following sequences:
1. Chromosome left portion (sequence 5' of TRC).
- The chromosome up to the TRC, ending with the sequence in components/left_linker.fa
2. Chromosome right portion (sequence 3' of TRC).
- The chromosome sequence after the TRC, starting with the sequence in components/right_linker.fa

Note:
- The contigs generated from the assembler likely include a portion of the TRC. Therefore, this should be removed.
- The following steps only account for assembling one contig, Contig1, into a singular chromosome. Repeat for each other contig.
   * If multiple contigs occur in the same chromosome (e.g. they are part of the same TRC), adapt the following workflow to concatenate them.
- The following steps assume that the TRC represented by Contig1 occurs in Chromosome 1 of the genome assembly. Replace Chromosome1 with any other number.

```bash
cp *path_to_left_chr_contig1_fa* genome/Contig1.chrleft.fa
cp *path_to_left_chr_contig1_fa* genome/Contig1.chrright.fa
```

For the fourth drafts generated in part 8, combine all of their annotations into one giant gff file, and simplify the chromosome ID. 
```bash
cat drafts/thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.gff drafts/thirddraft.V1.VAR.gff.thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.gff drafts/thirddraft.V1.VAR.firstdraft.gff.thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.gff drafts/thirddraft.V1.VAR.seconddraft.gff.thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.gff drafts/thirddraft.V1.features.gff.thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.gff | awk -v TRC_NAME=Contig1 'BEGIN {FS=OFS="\t"; print("##gff-version 3")} {if(substr($1,1,1) == "#") {next}; $1 = TRC_NAME; print}' | sort -k4,4n > genome/Contig1.gff
```

Assign segment numbers to each Variant annotation.
```bash
python3 09_labelSegments.py genome/Contig1.gff
```
Notes:
- This script assumes that feature annotations (originally from features.tandem.gff) are labeled with SEGMENT=1, SEGMENT=2 (etc) indicating which segment each feature belongs to.
- This script will also automatically add new Segment annotations that will provide the boundaries of segments. Then, it will label each Variant with which segment it belongs to. These are useful for downstream analysis in step 10.

Rename the fourth draft contig fasta file.
```bash
awk -v TRC_NAME=Contig1 'NR == 1 {print ">"TRC_NAME} NR >= 2 {print}' drafts/thirddraft.V1.short.Contig1.sam.VAR.tab.consensus.fa > genome/Contig1.fa
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 9b. Combine left portion with tandem repeat cluster.

```bash
python3 scripts/09_concatenateFasta.py genome/Contig1.chrleft.fa genome/Contig1.fa --annotations genome/Contig1.gff.labeled.gff
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 9c. Combine with right portion.

```bash
python3 scripts/09_concatenateFasta.py genome/Contig1.chrleft.fa.cat.Contig1.fa genome/Contig1.chrright.fa --dontupdatecoords --annotations genome/Contig1.gff.labeled.gff.cat.gff
```
Rename. Note: keep the CHR_NAME variable to exactly 4 characters.
```bash
awk -v CHR_NAME=Chr1 'NR == 1 {print ">"CHR_NAME} NR >= 2 {print}' genome/Contig1.chrleft.fa.cat.Contig1.fa.cat.Chromosome1right.fa > genome/Chromosome1.fa
awk -v CHR_NAME=Chr1 'BEGIN {FS=OFS="\t"; print("##gff-version 3")} {if(substr($1,1,1) == "#") {next}; $1 = CHR_NAME; print}'  genome/Contig1.gff.labeled.gff.cat.gff.cat.gff > genome/Chromosome1.draftgff.gff
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 9d. Map back reads to check junctions.

```bash
minimap2 -t *THREADS* --eqx -ax map-ont genome/Chromosome1.fa reads/long.fa.chop.oriented.fa > genome/Chromosome1.long.sam
minimap2 -t *THREADS* --eqx -ax map-hifi genome/Chromosome1.fa reads/short.fa.chop.oriented.fa > genome/Chromosome1.short.sam
```
Load the following files into visualization software:
- genome/Chromosome1.fa
- genome/Chromosome1.draftgff.gff
- genome/Chromosome1.short.sam
- genome/Chromosome1.long.sam

Check the read mappings for both short and long datasets where the Contigs connect to the rest of the chromosome.
- If there are any gaps in read mapping, go back to step 9a, and reprepare the chrleft and chrright sequences.
- If the read mappings cover the sequence correctly, proceed to 9e.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 9e. Finalization

Combine chromosomes with other chromosomes to complete the genome.\
Prepare a file with all genomic chromosomes that are correctly assembled from the original genome (i.e. containing no TRC).\
```bash
cp *path_to_correct_chromosomes_fa* genome/chromosomes.correct.fa
```
Replace \*...\* with a list of chromosomes assembled in step 9c, e.g. genome/Chromosome1.fa
```bash
cat genome/chromosomes.correct.fa *...* > genome/genome.fa
```

Map back *all reads* to the entire genome.
```
minimap2 -t *THREADS* --eqx -ax map-ont genome/genome.fa reads/long.fa > genome.long.sam
minimap2 -t *THREADS* --eqx -ax map-hifi genome/genome.fa reads/short.fa > genome.short.sam
```

Combining the annotations.\
Replace \*...\* with a list of chromosomes assembled in step 9c, e.g. genome/Chromosome1.draftgff.gff
```bash
cat *...* > genome/genome.tigerpaw.gff
```
**Manually edit the file genome/genome.tigerpaw.gff**. The Variant annotations will contain have a few innacuracies, for example:
- Some variants represent polishing of other variants. For example, retrotransposons will be listed as a very large insertion, with many other variants within that insertion. **Combine these into one large insertion**.
- In some rare cases, there may be multiple deletions at the same location. **Merge these variants**.
- In some rare cases, variants between drafts will cancel each other out (see below). **Delete these variants**.
```tsv
Chr1	Consensus	Variant	43558	43558	.	+	.	SEGMENT=5;TYPE="SNV";REF="G";NEW="A";
Chr1	Consensus	Variant	43558	43558	.	+	.	SEGMENT=5;TYPE="SNV";REF="A";NEW="G";
```

> This prepared genome contains novel sequence, and gives insight into the evolution of the tandem repeat cluster(s) within this organism. Annotation of this new genome with modern annotation tools may reveal new functions and/or contained genes within certain segment variants.


Optional step: ![WORKFLOW-10-ANALYSIS](WORKFLOW-10-ANALYSIS.md)
