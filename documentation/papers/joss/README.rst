:orphan: Avoid Sphinx warning about this file ...

JOSS paper
----------

This directory contains a file on rMarkdown_ format for submission of Ocellaris 2019.0.1 to the Journal of Open Source Software, http://joss.theoj.org/

To build the pdf you must have ``pandoc`` and ``pandoc-citeproc`` installed and then run::

    pandoc --filter pandoc-citeproc paper.md -o paper.pdf

You can also use a preconfigured Docker container that will use the same template as the Open Journals, see https://github.com/openbases/builder-pdf::

    docker run -v $PWD:/data openbases/openbases-pdf pdf

.. _rMarkdown: https://rmarkdown.rstudio.com/

Tormod Landet, February 2019
