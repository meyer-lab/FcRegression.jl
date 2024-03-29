
figures:
	mkdir -p ./output
	julia --thread 32 -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figureAll()'

figure%.svg:
	mkdir -p ./output
	julia --thread 32 -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figure$*()'

clean:
	rm -rf ./output/*.svg ./output/*.pdf
