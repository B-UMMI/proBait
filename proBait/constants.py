"""
Purpose
-------

This module contains variables used by proBait modules.

Code documentation
------------------
"""


summary_text = """# proBait report

The report has the following sections:
- **Configuration**: values passed to proBait\'s parameters.
- **Coverage statistics**: coverage statistics determined by mapping the final set of baits against each input.
- **Depth per position**: depth of coverage per position. Vertical dashed lines are contig boundaries and green markers along the x-axis are the start positions of baits that were generated to cover regions not covered by baits. Contigs are ordered based on decreasing length.
- **Depth values distribution**: distribution of depth of coverage values for each input (y-axis is shared with "Depth per position" plot in the same line).

If you have any question or wish to report an issue, please go to proBait\'s [GitHub](https://github.com/rfm-targa/proBait) repo.
"""
