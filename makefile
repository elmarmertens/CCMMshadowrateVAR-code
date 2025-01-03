
texcmd=pdflatex -interaction=batchmode

all  : CCMMpaperTablesAndFigures.pdf CCMMsupplementTablesAndFigures.pdf batch tablesandfigures

.PHONY : batch tablesandfigures clean

runbatch.m : preparebatch.sh
	sh preparebatch.sh goVAR*.m

batch : runbatch.m 
	sh goparbatch.sh runbatch.m

tablesandfigures : runtablesandfigures.m foo
	sh goseqbatch.sh runtablesandfigures.m

tex : CCMMpaperTablesAndFigures.pdf CCMMsupplementTablesAndFigures.pdf 

CCMMpaperTablesAndFigures.pdf : CCMMpaperTablesAndFigures.tex $(shell find foo -type f)
	$(texcmd) CCMMpaperTablesAndFigures.tex
	$(texcmd) CCMMpaperTablesAndFigures.tex
	$(texcmd) CCMMpaperTablesAndFigures.tex

CCMMsupplementTablesAndFigures.pdf : CCMMsupplementTablesAndFigures.tex $(shell find foo -type f)
	$(texcmd) CCMMsupplementTablesAndFigures.tex
	$(texcmd) CCMMsupplementTablesAndFigures.tex
	$(texcmd) CCMMsupplementTablesAndFigures.tex	

clean :
	rm -f CCMMpaperTablesAndFigures.pdf CCMMsupplementTablesAndFigures.pdf
	rm -f foo


