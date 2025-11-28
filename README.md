# proBait

proBait is a tool for designing baits for target enrichment/capture experiments. proBait starts by shredding a FASTA file with reference sequences to create an initial set of baits. This is followed by an iterative mapping approach (that uses [minimap2](https://github.com/lh3/minimap2)) in which the initial set of baits is mapped against the remaining input files to determine the sequence regions not covered by the baits and generate new baits to cover those regions according to the parameters defined by the user. proBait also includes options to cluster the generated bait set with [MMseqs2](https://github.com/soedinglab/MMseqs2) to reduce the redundancy of te bait set and the removal of baits that are similar to a host or contaminant by mapping the generated bait set against the host (e.g. human genome) and contaminant sequences (e.g. contaminants present in reagent kits).

### The [documentation](https://probait.readthedocs.io/en/latest/) includes detailed information about proBait's implementation, usage, and several tutorials.

## News

## Citation

Please cite this repository if you use proBait.

> Mamede, R., & Ramirez, M. (2024). proBait (Version 0.1.0) [Computer software]. https://github.com/B-UMMI/proBait
