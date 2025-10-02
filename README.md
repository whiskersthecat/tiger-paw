# tiger-paw

![tiger_paw logo](/assets/tiger_paw.png)



Workflow for very accurately closing and annotating gaps caused by large tandem repeat clusters (e.g. 400 copies of 10,000 base pair repeat) in eukaryotic genome assemblies.

Used to assemble 10 MB (10 million base pairs) and 4 MB of the two Nucleolus Organizer Regions (NORs) in [The Telomere to Telomere Gapless Lettuce (_lactuca sativa c. salinas_) Genome Assembly](https://kittishgames.com/pounce/).

**Background**: 
* In large tandem repeats, biological variation between consecutive repeat segments is minimal and less prevalent than sequencing errors, rendering the regions unresolveable by classical genome assembly algorithms employed by modern genome assemblers such as *Hifiasm* (String Overlap Graphs) or *Verkko* (DeBrujin Graphs).
* In many modern genome assemblies, these regions are left as gaps and filled with Ns or improperly scaffolded together.
* Correctly resolving these regions is very difficult, but gives tremendous insight into mechanisms of genomic evolution and greatly improves contiguity of the genome assembly.

**Notes**:
- This workflow is intended to produce a near-perfect sequence from an input of noisy long reads in combination with accurate shorter reads and a large time investment in manual alignment of **blocks** (basically what the tiger is doing).

