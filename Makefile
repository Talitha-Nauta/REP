#Comment if ploting normal or short 
all:
	latexmk -lualatex -shell-escape plot
	#latexmk -lualatex -shell-escape plot_short
	make -s clean

clean:
	rm -rf *.log *.out *.synctex.gz *.blg *.bbl *.auxlock *.nav *.snm *.toc *.fdb_latexmk *.fls *.aux