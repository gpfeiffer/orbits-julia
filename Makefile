# Makefile for "Computational Aspects of Complex Reflection Groups"
#
#   make             build the book (exports/book.pdf)
#   make clean       remove _book_build/ and exports/book.pdf
#
# nb2book.py executes the notebooks in memory (Julia kernel), moves the
# exercise solutions into an appendix, exports everything to LaTeX through
# templates/book, and compiles with xelatex.  The notebooks themselves are
# never modified.  See the docstring of nb2book.py for details.

# The notebooks in the book, as listed in nb2book.py (PREFACE, CHAPTERS).
# Drafts such as notebooks/computing.ipynb are deliberately left out.
NOTEBOOKS = $(addprefix notebooks/, $(addsuffix .ipynb, preface orbits coxeter enumerate linear))
TEMPLATE  = $(wildcard templates/book/*)
IMAGES    = $(wildcard images/*)

.PHONY: all clean

all: exports/book.pdf

exports/book.pdf: $(NOTEBOOKS) $(TEMPLATE) $(IMAGES) nb2book.py
	python3 nb2book.py

clean:
	rm -rf _book_build exports/book.pdf
