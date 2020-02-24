all: msAux ms

ms:
	pandoc -H style.sty -V fontsize=12pt -V geometry:margin=1in --bibliography=/Users/kylechezik/Documents/Reference_Literature/bibRef.bib --csl=apa.csl  ms.md -o ms.pdf --pdf-engine=pdflatex

msAux:
	pandoc -H style.sty -V fontsize=12pt -V geometry:margin=0.75in --bibliography=/Users/kylechezik/Documents/Reference_Literature/bibRef.bib --csl=apa.csl  ms_aux.md -o ms_aux.pdf

ms_docx:
	pandoc -H style.sty -V fontsize=12pt -V geometry:margin=1in --bibliography=/Users/kylechezik/Documents/Reference_Literature/bibRef.bib --csl=apa.csl  ms.md -o ms.docx

