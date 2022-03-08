
figures:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figureAll()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figure$*()'

clean:
	rm -rf *.svg *.pdf
