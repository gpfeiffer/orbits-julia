# Makefile for "Computational Aspects of Complex Reflection Groups"
#
#   make             build the book (exports/book.pdf)
#   make clean       remove _book_build/ and exports/book.pdf
#
# nb2book.py executes the notebooks in memory (Julia kernel), moves the
# exercise solutions into an appendix, exports everything to LaTeX through
# templates/book, and compiles with xelatex.  The notebooks themselves are
# never modified.  See the docstring of nb2book.py for details.

NOTEBOOKS = $(wildcard notebooks/*.ipynb)
TEMPLATE  = $(wildcard templates/book/*)

.PHONY: all clean

all: exports/book.pdf

exports/book.pdf: $(NOTEBOOKS) $(TEMPLATE) nb2book.py
	python3 nb2book.py

clean:
	rm -rf _book_build exports/book.pdf
