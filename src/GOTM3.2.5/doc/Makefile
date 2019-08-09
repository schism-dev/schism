#$Id: Makefile,v 1.10 2005-08-11 13:59:25 lars Exp $
#
# Makefile for making the source code documentation. 
#

SHELL	   = /bin/sh

SRCDIR    = ../src

PROTEXSRC  = gotm.tex              \
             meanflow.tex          \
             turbulence.tex        \
             airsea.tex            \
             observations.tex      \
             output.tex            \
             util.tex              \
             seagrass.tex          \
#             sediment.tex          \
#             bio.tex

TEXSRC     = introduction.tex      \
             mainIntro.tex         \
             meanflowIntro.tex     \
             turbulenceIntro.tex   \
             airseaIntro.tex       \
             observationsIntro.tex \
             outputIntro.tex       \
             utilIntro.tex         \
             extraIntro.tex        \
             cases.tex             \
             gotm.bib


TEXHTMLSRC = introduction.tex      \
             mainIntro.tex         \
             meanflowIntroHTML.tex \
             turbulenceIntro.tex   \
             airseaIntro.tex       \
             observationsIntro.tex \
             outputIntro.tex       \
             utilIntroHTML.tex     \
             extraIntro.tex        \
             cases.tex             \
             gotm.bib



.PHONY :     all clean realclean distclean protex \
             a4_dvi letter_dvi html_dvi           \
             a4_ps  letter_ps  html_ps            \
             a4_pdf letter_pdf html_pdf           \
             html_web


all:         a4_pdf letter_pdf html_pdf html_web


# do the pdf-versions
a4_pdf:      protex a4.pdf

letter_pdf:  protex letter.pdf

html_pdf:    protex html.pdf

# do the postscript-versions
a4_ps:       protex a4.ps

letter_ps:   protex letter.ps

html_ps:     protex html.ps

# do the dvi-versions
a4_dvi:      protex a4.dvi

letter_dvi:  protex letter.dvi

html_dvi:    protex html.dvi

# do the web-version
html_web:    protex html

protex:
	$(MAKE) realclean -C $(SRCDIR)
	$(MAKE) doc       -C $(SRCDIR)

a4.dvi: a4.tex  $(TEXSRC) $(PROTEXSRC) figures
	latex  a4.tex
	bibtex a4
	latex  a4.tex
	latex  a4.tex
	latex  a4.tex

letter.dvi: letter.tex $(TEXSRC) $(PROTEXSRC) figures
	latex  letter.tex
	bibtex letter
	latex  letter.tex
	latex  letter.tex
	latex  letter.tex

html.dvi: html.tex $(TEXSRC) $(PROTEXSRC) figures
	latex  html.tex
	bibtex html
	latex  html.tex
	latex  html.tex
	latex  html.tex

a4.ps: a4.dvi
	dvips -o $@ a4.dvi

letter.ps: letter.dvi
	dvips -o $@ letter.dvi

html.ps: html.dvi
	dvips -o $@ html.dvi

a4.pdf: a4.ps
	ps2pdf a4.ps

letter.pdf: letter.ps
	ps2pdf letter.ps

html.pdf: html.ps
	ps2pdf html.ps

html:   html.dvi
	latex2html -antialias_text -white html.tex
	touch html

pdf_scp: a4.pdf letter.pdf
	scp -r a4.pdf letter.pdf gotm@80.198.12.199:doc/
	touch pdf_scp

html_scp: html
	scp -r html/* gotm@80.198.12.199:html/
	touch html_scp


figures:
	wget -nd -nH -q http://www.gotm.net/data_files/doc_figures.tar.gz 
	tar -xvzf doc_figures.tar.gz
	$(RM) doc_figures.tar.gz

clean:
	$(RM) *.log *.toc *.aux *.bbl *.blg *.bac *~

realclean: clean
	$(RM) *.dvi *.ps *.pdf
	$(RM) -rf html/ 
	$(RM) -rf $(PROTEXSRC)

distclean: realclean
	$(RM) *.gz *.zip *.Z
	$(RM) html html_scp pdf_scp 
	$(RM) -rf figures

