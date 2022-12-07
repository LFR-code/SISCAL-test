export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
pdf:  AtlHerring_4XBOF.pdf AtlHerr_TandF.pdf supplemental.pdf
html:  AtlHerring_4XBOF.html AtlHerr_TandF.html supplemental.html
docx:  AtlHerring_4XBOF.docx supplemental.docx

all: pdf html docx


## pdf
%.pdf: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", clean = TRUE, output_format = "bookdown::pdf_document2")'

# html
%.html: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", clean =TRUE, output_format = "bookdown::html_document2")'

# docx
%.docx: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", clean = TRUE, output_format = "word_document")'

# all formats (depends on yaml header)
%.all: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", clean = TRUE, output_format = "all")'


## appendix.tex
# appendix.tex: appendix.md makefile
# 	$(PANDOC) -s $(NOPRETMP) $(ABV) $(BIB) $(BIBPRE) $(CSL) $(FONT) $(TIMS) -o $@ $<

# TandF
# TandF.pdf: TandF.tex makefile
# 	$(LATEX) $<

# TandF.tex: TandF.md makefile
# 	$(PANDOC) -s -S $(PTMP) -o $@ $<	

# TandFnopre.tex: TandF.md makefile
# 	$(PANDOC) -s -S $(NOPRETMP) -o $@ $<	

## cleanup
clean:
	rm -f *.{pdf,doc,docx,tex,aux,log,out,toc,bbl}
