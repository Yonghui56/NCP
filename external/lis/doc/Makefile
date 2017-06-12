# Makefile for LaTeX Systems

DOC=lis-manual

LATEX=latex
JLATEX=platex

DVIPDF=dvipdfmx

all: en ja 

en: $(DOC)-en.pdf

ja: $(DOC)-ja.pdf

%.pdf: %.dvi
	$(DVIPDF) $< 

%-en.dvi: %-en.tex
	$(LATEX) $<
	$(LATEX) $<

%-ja.dvi: %-ja.tex
	$(JLATEX) $<
	$(JLATEX) $<

clean: 
	-$(RM) *.dvi *.aux *.toc *.log

distclean: clean
	-$(RM) *.pdf
