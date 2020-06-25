
all: figures depletion temporal output/depletion/manuscript.html output/translation/manuscript.html

venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uqr requirements.txt
	touch venv/bin/activate

depletion: figure2.svg figure3.svg

figures:
	julia -e 'using Pkg; Pkg.activate("."); using FcgR; FcgR.figureAll()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using FcgR; FcgR.figure$*()'

temporal:
	julia -e 'using Pkg; Pkg.activate("."); using FcgR; FcgR.plotTemporal()'

output/%/manuscript.md: venv manuscripts/%/*.md
	mkdir -p ./output/$*
	. venv/bin/activate && manubot process --content-directory=manuscripts/$*/ --output-directory=output/$*/ --log-level=WARNING

output/%/manuscript.html: venv output/%/manuscript.md figureB1.svg figureB2.svg figureB3.svg
	cp *.svg output/$*/
	. venv/bin/activate && pandoc \
		--from=markdown --to=html5 --filter=pandoc-fignos --filter=pandoc-eqnos --filter=pandoc-tablenos \
		--bibliography=output/$*/references.json \
		--csl=common/templates/manubot/style.csl \
		--metadata link-citations=true \
		--include-after-body=common/templates/manubot/default.html \
		--include-after-body=common/templates/manubot/plugins/table-scroll.html \
		--include-after-body=common/templates/manubot/plugins/anchors.html \
		--include-after-body=common/templates/manubot/plugins/accordion.html \
		--include-after-body=common/templates/manubot/plugins/tooltips.html \
		--include-after-body=common/templates/manubot/plugins/jump-to-first.html \
		--include-after-body=common/templates/manubot/plugins/link-highlight.html \
		--include-after-body=common/templates/manubot/plugins/table-of-contents.html \
		--include-after-body=common/templates/manubot/plugins/lightbox.html \
		--mathjax --variable math="" \
		--include-after-body=common/templates/manubot/plugins/math.html \
		--include-after-body=common/templates/manubot/plugins/hypothesis.html \
		--output=output/$*/manuscript.html output/$*/manuscript.md

coverage.cob:
	julia -e 'using Pkg; Pkg.add("Coverage"); using Coverage; Pkg.activate("."); Pkg.test("FcgR"; coverage=true); coverage = process_folder(); LCOV.writefile("coverage-lcov.info", coverage)'
	pip3 install --user lcov_cobertura
	python3 ~/.local/lib/python3.8/site-packages/lcov_cobertura.py coverage-lcov.info -o coverage.cob

clean:
	rm -rf *.svg venv output
