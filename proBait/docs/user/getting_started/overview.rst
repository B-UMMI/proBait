Overview
========

About proBait
-------------

proBait is a tool for designing baits for target enrichment/capture experiments. An overview of 
proBait's workflow is shown below.

.. image:: /_static/images/proBait_overview.png
   :width: 1000px
   :align: center

proBait starts by breaking a FASTA file with reference sequences (if a file containing the paths to the reference 
FASTA files is provided through the `--refs` option, otherwise it uses the first FASTA file as reference) to create 
an initial set of baits. For proBait to generate baits, it is necessary to use the `--generate-baits` option, 
otherwise proBait will only perform a bait coverage evaluation based on a set of baits provided by the user through
the `--baits` option. After generating the initial set of baits based on the reference sequences, proBait use an 
iterative mapping approach (using `minimap2 <https://github.com/lh3/minimap2>`_) in which the initial set of baits 
is mapped against the remaining input files to determine the sequence regions not covered by the baits and generate 
new baits to cover those regions according to the parameters defined by the user. proBait also includes options to 
cluster the generated bait set with `MMseqs2 <https://github.com/soedinglab/MMseqs2>`_ to reduce the redundancy of 
the bait set and the removal of baits that are similar to a host or contaminant by mapping the baits against the host
(e.g. human genome) and contaminant sequences (e.g. contaminants present in reagent kits).

The option to provide a previosuly defined set of baits allows users to evaluate the coverage of their baits against
a set of target sequences, which can be useful to assess the performance of existing bait sets (generated with proBait 
or other methods) or to compare different bait sets for a given set of target sequences.

The coverage evaluation can be performed for multiple identity and coverage (`-ri` and `-rc` options, respectively) 
values to assess the performance of the baits under different parameter combinations.

Citation
--------

When using proBait please cite its GitHub repository:

```
Mamede, R., & Ramirez, M. (2024). proBait (Version 0.1.0) [Computer software]. https://github.com/B-UMMI/proBait
```

Licensing
---------

This project is licensed under the `GPLv3 license <https://github.com/B-UMMI/proBait/blob/main/LICENSE>`_.
The source code of proBait is available on `GitHub <https://github.com/B-UMMI/proBait>`_.

Funding
-------
