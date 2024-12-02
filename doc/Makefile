#  set the lyx location here, if it is not in your path:
#LYX=/Applications/LyX.app/Contents/MacOS/lyx
LYX=lyx

SHELL = /bin/sh
.SUFFIXES: .txt .ps .html .lyx .pdf .tex .pic

.lyx.pdf:
	$(LYX) --export pdf2 $<
		
.lyx.tex:
	$(LYX) --export pdflatex $<
		
.lyx.txt:
	$(LYX) --export text $<
		
.lyx.ps:
	$(LYX) --export ps $<

.pic.pdf: $<
	groff -p $< | ps2pdf - $@

AmberClassic.pdf: AmberClassic.lyx amberclassic.bib
	-$(LYX) --export pdf2 AmberClassic.lyx

clean::
	/bin/rm -f AmberClassic.pdf
