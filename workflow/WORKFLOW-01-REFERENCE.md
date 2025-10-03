# 1. Create Reference

Overview: generate a highly quality generic repeat segment reference, and a set of feature annotations.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 1a. Generate reference segment sequence.

```bash
mkdir reference
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

### 1a. Generate reference segment feature annotations.

Annotations are useful in visualizing locations of variants.\
```
echo "##gff-version 3" reference/features.gff
```
Blast the reference sequence against the NCBI database:
```tsv
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
```
Sort the hits by the highest percent identity.\

Extract useful hits and manually write them into the *features.gff* file.\
For example, this annotations file was used in the assembly of the NOR region in lettuce:
```tsv
##gff-version 3
Ls_NOR_Ref	annotation	rRNA	1	1574	.	+	.	TYPE=Ls_28S;
Ls_NOR_Ref	annotation	rRNA	5681	7490	.	+	.	TYPE=Ls_18S;
Ls_NOR_Ref	annotation	rRNA	7747	7902	.	+	.	TYPE=Ls_5S;
Ls_NOR_Ref	annotation	rRNA	8132	9847	.	+	.	TYPE=Ls_28S;
Ls_NOR_Ref	annotation	spacer	2052	3200	.	+	.	TYPE=spacer;
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-02-READS](WORKFLOW-02-READS.md)