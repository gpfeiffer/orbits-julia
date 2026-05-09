.SUFFIXES: .pdf .png


.tex.pdf:
	pdflatex $<

.pdf.png:
	convert -density 300 $< $@

all: e7tree.png e7direct.png bfsdfs.png orbit.png

clean:
	rm -f *.aux *.log

distclean: clean
	rm -f *.png


