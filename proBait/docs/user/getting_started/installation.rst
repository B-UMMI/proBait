Installation
============

Conda
.....

Install the latest released version using :

::

	conda create -c bioconda -c conda-forge -n proBait ""

If you're having issues installing proBait through conda, please verify that you are using
conda>=22.11, and enable the libmamba solver, which might speed up the installation process.
You can also install `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_ and run the following command:

::

	mamba create -c bioconda -c conda-forge -n proBait ""

Pip
...

Install using `pip <https://pypi.org/project/proBait/>`_:

::

	pip3 install proBait


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
