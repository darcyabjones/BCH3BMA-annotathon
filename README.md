BCH3BMA-annotathon
==================

Project Makefiles and LaTeX/knitr documents for a university meta-genomics analysis assignment.

I am a stubborn budding biologist/bioinformatician and decided to make my life significantly harder by attempting to make my homework reproducible. Most analysis can be done entirely using the makefile commands.mk. Swissprot structural and functional analysis was done using EMBL and ExPASy webserver tools. Phylogenetic tree prettifying also has to be done in a GUI interface such as dendroscope or figtree.

All files required for generating the main information files, including blast results, the phylogeny and sweave document can be made by running make in the working directory, ie:

	$ make

The Cladogram and VMD structure figures could not be generated exclusively from the command line and so they are included as pdfs.

Required software
-----------------

Emboss command line tools

Blast+

Entrez-direct

python3 with optional package biopython installed

T-Coffee

RAxML 8

TeXLive with optional package TeXshade installed

R with optional packages xtable and knitr installed
