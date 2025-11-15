![tiger_paw logo](assets/tiger_paw.png)

## Description

Workflow for very accurate closing, variant annotation, and analysis of gaps in eukaryotic genome assemblies caused by very large tandem repeat clusters (e.g. hundreds of back-to-back duplications of a ~10,000 base pair gene).

Used to assemble 10 MB (10 million base pairs) and 4 MB of the two Nucleolus Organizer Regions (ribosomal rRNA gene clusters) in [The Telomere to Telomere Gapless Lettuce (_lactuca sativa c. salinas_) Genome Assembly](https://kittishgames.com/pounce/).

![example_assembly](assets/example_assembly.png)

## Background

* In large tandem repeat clusters, subtle biological sequence variants between consecutive repeat segments are obscured by sequencing errors, rendering consecutive segments virtually indistinguishable.
* Modern genome assemblers like *Hifiasm* (string overlap graphs) or *Verkko* (Debrujin Graphs) cannot resolve these regions and leave gaps in the genome assembly.
* Correctly resolving these regions is very difficult, but greatly improves contiguity of the genome assembly and gives tremendous insight into mechanisms of genome evolution.

## Notes

- This workflow is intended to produce a near-perfect sequence from an input of noisy long reads in combination with accurate shorter reads and a large time investment in manual stacking of **blocks** (basically what the tiger is doing).

## Usage

The ***tiger-paw*** workflow is seperated into distinct steps (see below). Each step has an associated markdown file which describes exactly what scripts to run and analysis to perform.

Clone this repository.
```bash
git clone whiskersthecat/tiger-paw
```

Start with [WORKFLOW-00-INTRODUCTION](workflow/WORKFLOW-00-INTRODUCTION.md), which describes how to install the few software dependencies and the required input data.


## Overview
![tiger_paw workflow](assets/tigerpaw_workflow.png)
