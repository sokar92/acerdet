OUT= AcerDET-2.0

clean:
	rm $(OUT).aux
	rm $(OUT).dvi
	rm $(OUT).log

tex:                                                     
	latex $(OUT).tex 
	dvips $(OUT).dvi -o 
	gv $(OUT).ps 
