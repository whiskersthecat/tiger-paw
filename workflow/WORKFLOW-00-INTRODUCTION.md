# 0. Introduction


These **definitions** are used throughout the workflow.

* ***segment***: one repeat copy of a tandem repeat. A single long read should span many segments.
* ***chimera***: a single read in which two different (unrelated) molecules are concatenated together.
* ***duplex***: a single read in which the same molecule is read twice, with forward and reverse strand concatenated together.
* ***TRC***: Tandem Repeat CLuster.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


**tiger-paw** requires two large DNA sequencing datasets, ideally both with coverage at least 50x. For *long.fa*, the more the better.

* ***short.fa***. Reads are up to 30,000 base pairs, each read spans at least one segment, and sequence is **very accurate**. Reads are assumed to contain no chimeras or duplexes. For example, *Pacbio HiFi Reads*. 
> These reads will be used to generate common Haplotypes that will be used in assembly.

* ***long.fa***. Reads are longer than 30,000 base pairs and can span at **least three segments, ideally up to 20 or more**. They are assumed to be innacurate, and contain duplexes and potentially  chimeras. For example, *Oxford Nanopore Reads*. 
> These reads will be used to actually assemble the tandem repeat cluster. 

Note: in the case there are multiple datasets of either type, it is best to keep them seperate and run them in parallel. This helps confirm phenomena that is replicated among all datasets.

**tiger-paw** also requires an input of a reference segment. The reference segment does not have to be perfect, but ideally a similar length as a single segment and free of errors.




<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Dependencies:
- minimap2 (https://github.com/lh3/minimap2)
- blast Version 2.16.0 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)
   - Note: different blast Versions have very different performance. Newest blast version 2.17.0 struggled in identifying duplexes in HiFi reads.
- bedtools
- clustalw version 2.1
- python version 3.11

You can install the dependencies using mamba:
```bash
mamba create -n tigerpaw -c bioconda minimap2 blast=2.16.0 bedtools python=3.11 clustalw=2.1
```

Supplementary Software:\
It is important to use one of the following visualization software to visualize alignments on multiple steps:
- Interactive Genomics Viewer (https://igv.org/)
- CLC Genomics Workbench (https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/qiagen-clc-genomics/?cmpid=QDI_GA_DISC_CLC_SA&gad_source=1&gad_campaignid=6942460474&gbraid=0AAAAADbyWl0yf7q1l0zt_ZecsD7avy_bK&gclid=Cj0KCQjwovPGBhDxARIsAFhgkwSB0_xW4BRg1U7Z2qqdnT1Qtfd1Hrt5fg4xcGbaECuH1LFYdof3tpEaAsMCEALw_wcB)
- Any other modern genomics browser or visualization software.

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>


### 0A. Copy reads

Note: in all commands presented in the workflow, text in \* \* should be replaced with user parameters.

```bash
mkdir reads
cp *path_to_short_reads* reads/short.fa
cp *path_to_long_reads* reads/long.fa
```

<sub>\--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- </sub>

Next step: ![WORKFLOW-01-REFERENCE](WORKFLOW-01-REFERENCE.md)
