OUTNAME = compression-$(shell date +%F).pdf

all: compression.pdf

compression.pdf: q-compression.tex
	pdflatex q-compression.tex; bibtex q-compression; pdflatex q-compression.tex; pdflatex q-compression.tex;
	mv q-compression.pdf ${OUTNAME}

clean:
	rm -rf q-compression*.pdf
	rm -f q-compression.log 
	rm -f q-compression.dvi 
	rm -f q-compression.aux 
	rm -f q-compression.bbl
	rm -f q-compression.blg
	rm -f q-compression.bcf