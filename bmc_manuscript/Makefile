OUTNAME = compression-$(shell date +%F).pdf

all: compression.pdf

compression.pdf: manuscript.tex
	pdflatex manuscript.tex; bibtex manuscript; pdflatex manuscript.tex; pdflatex manuscript.tex;
	mv manuscript.pdf ${OUTNAME}

clean:
	rm -rf q-compression*.pdf
	rm -f q-compression.log
	rm -f q-compression.dvi
	rm -f q-compression.aux
	rm -f q-compression.bbl
	rm -f q-compression.blg
	rm -f q-compression.bcf
