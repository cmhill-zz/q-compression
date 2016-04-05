OUTNAME = compression-$(shell date +%F).pdf

all: compression.pdf

compression.pdf: bmc_article.tex
	pdflatex bmc_article.tex; bibtex bmc_article; pdflatex bmc_article.tex; pdflatex bmc_article.tex;
	mv bmc_article.pdf ${OUTNAME}

clean:
	rm -rf q-compression*.pdf
	rm -f q-compression.log
	rm -f q-compression.dvi
	rm -f q-compression.aux
	rm -f q-compression.bbl
	rm -f q-compression.blg
	rm -f q-compression.bcf
