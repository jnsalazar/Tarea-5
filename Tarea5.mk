Results_hw5.pdf: curva_ajuste.pdf
	pdflatex Results_hw5.tex
curva_ajuste.pdf: fit.txt
	python Plots.py
fit.txt: a.out
	./a.out
a.out: RadialVelocities.dat
	cc CurvaRotacion.c