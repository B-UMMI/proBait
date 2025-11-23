# proBait

proBait is a tool for designing baits/probes for target enrichment/target capture experiments of bacterial pathogens. proBait starts by shredding a set of reference genome assemblies (in FASTA format) to create an initial set of baits. This is followed by an iterative mapping approach (that uses [minimap2](https://github.com/lh3/minimap2)) in which the initial set of baits is mapped against the input genome assemblies to determine the genomic regions not covered by the baits and generate new baits to cover those regions according to the parameters set by the user. proBait also includes options to cluster the generated bait set with [MMseqs2](https://github.com/soedinglab/MMseqs2) to reduce redundancy and the removal of baits that are similar to a host or contaminant by mapping the generated bait set against the host and contaminant sequences.

## Installation

We recommend creating a separate environment with conda to install proBait and its dependencies. You can do this by installing [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install).

### Pip

```
pip3 install probait
```

### Source

Just `cd` into the `proBait` directory after cloning this repository and run:

```
pip install .
```

#### Python dependencies

- biopython>=1.83
- plotly>=5.22.0
- pandas>=1.5.3
- datapane>=0.17.0

These dependencies are automatically installed in any of the aforementioned installation methods.

#### Other dependencies

- [minimap2](https://github.com/lh3/minimap2)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)

These dependencies are not installed automatically. Please install them in the environment you are working on.

## Usage

```
usage: proBait.py [-h] -i INPUT_FILES -o OUTPUT_DIRECTORY [-gb] [-b BAITS] [-bp BAIT_PROPORTION] [-rf REFS] [-msl MINIMUM_SEQUENCE_LENGTH]
                  [-bs BAIT_SIZE] [-bo BAIT_OFFSET] [-bi BAIT_IDENTITY] [-bc BAIT_COVERAGE] [-mr MINIMUM_REGION] [-me MINIMUM_EXACT_MATCH]
                  [-c] [-ci CLUSTER_IDENTITY] [-cc CLUSTER_COVERAGE] [-e EXCLUDE] [-ep EXCLUDE_PIDENT] [-ec EXCLUDE_COVERAGE] [-t THREADS]
                  [-r] [-ri REPORT_IDENTITIES [REPORT_IDENTITIES ...]] [-rc REPORT_COVERAGES [REPORT_COVERAGES ...]] [-tsv]

Purpose
-------

Generate baits for target capture experiments.

Code documentation
------------------

options:
  -h, --help            show this help message and exit
  -i INPUT_FILES, --input-files INPUT_FILES
                        Path to the directory that contains the input FASTA files.
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Path to the output directory that will be created to store output files. The process will exit if the directory
                        already exists.
  -gb, --generate-baits
                        Pass this parameter to generate baits based on the sequences in the input files.
  -b BAITS, --baits BAITS
                        Path to a FASTA file with baits generated on a previous run.
  -bp BAIT_PROPORTION, --bait-proportion BAIT_PROPORTION
                        Path to a TSV file with data about bait proportion to use when evaluating bait performance and generating the
                        reports.
  -rf REFS, --refs REFS
                        Path to a file with the basenames of the input FASTA files, one basename per line, that will be used as references
                        to create the initial set of baits. The references are shredded into baits according to the --bait-size and
                        --bait-offset values.
  -msl MINIMUM_SEQUENCE_LENGTH, --minimum-sequence-length MINIMUM_SEQUENCE_LENGTH
                        Do not generate baits for sequences shorter than this value.
  -bs BAIT_SIZE, --bait-size BAIT_SIZE
                        The length of the baits in bases. All the baits that are generated during the process have a size equal to the
                        passed value.
  -bo BAIT_OFFSET, --bait-offset BAIT_OFFSET
                        Start position offset between consecutive baits.
  -bi BAIT_IDENTITY, --bait-identity BAIT_IDENTITY
                        Minimum percent identity to accept an alignment between a bait and a region of an input sequence.
  -bc BAIT_COVERAGE, --bait-coverage BAIT_COVERAGE
                        Minimum proportion of a bait that must align to accept an alignment.
  -mr MINIMUM_REGION, --minimum-region MINIMUM_REGION
                        The process will only generate new baits for uncovered regions with length greater than this value.
  -me MINIMUM_EXACT_MATCH, --minimum-exact-match MINIMUM_EXACT_MATCH
                        Minimum number of N sequential matching bases in an alignment to accept it.
  -c, --cluster         Cluster set of baits to remove similar baits and reduce redundancy.
  -ci CLUSTER_IDENTITY, --cluster-identity CLUSTER_IDENTITY
                        Exclude baits with an identity value to the cluster representative equal to or higher than this value.
  -cc CLUSTER_COVERAGE, --cluster-coverage CLUSTER_COVERAGE
                        Exclude baits with a coverage value to the cluster representative equal to or higher than this value.
  -e EXCLUDE, --exclude EXCLUDE
                        Path to a FASTA file containing sequences to which baits must not be specific.
  -ep EXCLUDE_PIDENT, --exclude-pident EXCLUDE_PIDENT
                        Exclude baits with an identity value to a region of a sequence to exclude equal or higher than this value.
  -ec EXCLUDE_COVERAGE, --exclude-coverage EXCLUDE_COVERAGE
                        Exclude baits with a coverage value to a region of a sequence to exclude equal or higher than this value.
  -t THREADS, --threads THREADS
                        Number of threads passed to minimap2 and MMseqs2.
  -r, --report          Evaluate bait performance against input sequences and generate an interactive report with results.
  -ri REPORT_IDENTITIES [REPORT_IDENTITIES ...], --report-identities REPORT_IDENTITIES [REPORT_IDENTITIES ...]
                        List of identity values used to evaluate bait performance. proBait will generate a report per identity value. An
                        equal number of coverage values must be provided to the --report-coverages parameter to pair with the identity
                        values.
  -rc REPORT_COVERAGES [REPORT_COVERAGES ...], --report-coverages REPORT_COVERAGES [REPORT_COVERAGES ...]
                        List of coverage values used to evaluate bait performance. proBait will generate a report per coverage value. An
                        equal number of identity values must be provided to the --report-identities parameter to pair with the coverage
                        values.
  -tsv, --tsv-output    Output bait set in TSV format (the first column includes the bait sequence identifier, and the second column includes
                        the bait DNA sequence).
```

## Citation

Please cite this repository if you use proBait.

> Mamede, R., & Ramirez, M. (2024). proBait (Version 0.1.0) [Computer software]. https://github.com/B-UMMI/proBait
