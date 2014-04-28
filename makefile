# Darcy's standard makefile for Sweave documents. 

INFILE = projectc
SEQ=${PWD}/sequences
BLAST=${PWD}/blast
PHYLO=${PWD}/phylo
FUNCT=${PWD}/function

# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.

.PHONY: $(INFILE).pdf all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(INFILE).pdf

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# you might have.

#$(INFILE).Rnw: $(INFILE)_blast_results.xml $(INFILE)_blast_results_summary.bib $(INFILE)_blast_results_summary.txt

${INFILE}.Rnw: ${SEQ}/sequence_longest_orf.fa ${BLAST}/sequence_longest_orf.tsv ${PHYLO}/best_tree.nwk ${BLAST}/sequences_msa.fa ${FUNCT}/sequence_longest_orf_cdd.tsv ${FUNCT}/sequence_longest_orf_swissprot.tsv


# Doing references like this incase i want to include some references via straight bibtex.
${INFILE}.bib: ~/tex/bibtex/annotathon.bib
	cp $< $@

%.Rnw: %.bib

include commands.mk

%.tex: %.Rnw
		Rscript -e "library(knitr); knit(\"$<\")"

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interactive=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(INFILE).pdf: $(INFILE).tex
		latexmk -pdf -pdflatex="pdflatex -interactive=nonstopmode" -use-make $(INFILE).tex
		latexmk -c
