VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
#VERSION = git
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)
LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  all
all:  defs # html pdf
pdf: $(PDFS)
html: user_guide

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/
	- rm -rf input_xx.xsl
	- rm -f *.html *.txt *.xml
	- rm -rf ../../Doc/INPUT_ALL_CURRENTS.*


user_guide: user_guide.pdf
	- rm -rf user_guide/
	latex2html \
                -t "User's Guide for the energy_current package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                user_guide.tex
	cd user_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in user_guide/user_guide.html"
	@echo "***"
	@echo ""


defs: input_xx.xsl INPUT_ALL_CURRENTS.html INPUT_ALL_CURRENTS.txt link_on_main_doc

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1; \
	fi) ; fi)

INPUT_ALL_CURRENTS.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_ALL_CURRENTS.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

link_on_main_doc:
	-@( cd ../../Doc ; ln -fs ../QEHeat/Doc/INPUT_ALL_CURRENTS.html . ; \
	ln -fs ../QEHEAT/Doc/INPUT_ALL_CURRENTS.xml . ; \
	ln -fs ../QEHeat/Doc/INPUT_ALL_CURRENTS.txt .) 
