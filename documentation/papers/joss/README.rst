:orphan: Avoid Sphinx warning about this file ...

JOSS paper
----------

This directory contains a file on rMarkdown_ format for submission of Ocellaris 2019.0.1 to the Journal of Open Source Software, http://joss.theoj.org/

To build the pdf you must have ``pandoc`` and ``pandoc-citeproc`` installed and then run::

    pandoc --filter pandoc-citeproc paper.md -o paper.pdf

.. _rMarkdown: https://rmarkdown.rstudio.com/

Tormod Landet, February 2019
