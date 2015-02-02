OUTNAME = compression-$(shell date +%F).pdf

all: compression.pdf

compression.pdf: document.tex
	pdflatex document.tex; bibtex document; pdflatex document.tex; pdflatex document.tex;
	mv document.pdf ${OUTNAME}

clean:
	rm -rf document*.pdf
	rm -f document.log 
	rm -f document.dvi 
	rm -f document.aux 
	rm -f document.bbl
	rm -f document.blg
	rm -f document.bcf