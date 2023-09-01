# Makefile
DOC_FILES = $(shell find ../../.. -type f \( -iname \*.md -o -iname \*.bib \))
GENFOAM_SRC = $(shell find ../../../GeN-Foam -type f \( -iname \*.H -o -iname \*.C \))


.PHONY: all
all: html

.PHONY: html
html: html/index.htm

html/index.htm : Doxyfile $(DOC_FILES) $(GENFOAM_SRC) 
	@doxygen
	@touch $@

all :
	@cd ../../ ; ln -s doxygen/Doxygen/html/index.html Documentation_doxygen

clean :
	@echo Deleting existing documentation
	@cd ../../ ; rm -v Documentation_doxygen
	@rm -rv html DTAGS 

