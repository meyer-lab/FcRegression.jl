
figures:
	julia --thread 32 -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figureAll()'

figure%.svg:
	julia --thread 32 -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figure$*()'

clean:
	rm -rf *.svg *.pdf
