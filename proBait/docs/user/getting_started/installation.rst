Installation
============

Conda
.....

Install the latest released version using :

::

	conda create -c bioconda -c conda-forge -n proBait "python=3.11"

Followed by the following commands to install proBait and its dependencies:

::

	conda activate proBait
	conda install -c bioconda -c conda-forge biopython plotly pandas datapane minimap2 mmseqs2
	pip3 install probait

If you're having issues installing proBait through conda, please verify that you are using
conda>=22.11, and enable the libmamba solver, which may speed up the installation process.

Python dependencies
...................

* biopython >=1.83
* plotly >=5.22.0
* pandas >=1.5.3
* datapane >=0.17.0

.. note::
	These dependencies are defined in the `requirements <https://github.com/B-UMMI/proBait/blob/main/proBait/requirements.txt>`_
	file and should be automatically installed when using conda or pip.

Other dependencies
..................

* `minimap2 <https://github.com/lh3/minimap2>`_
* `MMseqs2 <https://github.com/soedinglab/MMseqs2>`_

.. important::
	Installation through conda should take care of all dependencies. If you install through
	pip you need to ensure that you have minimap2 and MMseqs2 installed.
