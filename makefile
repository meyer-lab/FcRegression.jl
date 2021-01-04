
all: figures output/manuscript.html

venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uqr requirements.txt
	touch venv/bin/activate

figures:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figureAll()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using FcRegression; FcRegression.figure$*()'

output/manuscript.md: venv manuscripts/*.md
	mkdir -p ./output
	. venv/bin/activate && manubot process --content-directory=manuscripts --output-directory=output --cache-directory=cache --skip-citations --log-level=INFO
	git remote rm rootstock

output/manuscript.html: venv output/manuscript.md
	cp *.svg output/
	. venv/bin/activate && pandoc --verbose \
		--defaults=./common/templates/manubot/pandoc/common.yaml \
		--defaults=./common/templates/manubot/pandoc/html.yaml output/manuscript.md

clean:
	rm -rf *.svg venv output
