
LATEX=latex
DVIPS=dvips

data_DATA=autotools_flow.eps fem_example.eps

.tex.eps:
	${LATEX} $< && ${DVIPS} -E -o $@ $*.dvi && rm -f $*.dvi

mostlyclean-local:
	rm -f *.dvi *.eps

clean-local:
	rm -f *.log *.aux
