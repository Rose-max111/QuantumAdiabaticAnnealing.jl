JL = julia --project

default: init test

init:
	$(JL) -e 'using Pkg; Pkg.precompile(); Pkg.activate("project"); Pkg.develop(path="."), Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile(); Pkg.activate("project"); Pkg.update(); Pkg.precompile()'

test:
	$(JL) -e 'using Pkg; Pkg.test()'

coverage:
	$(JL) -e 'using Pkg; Pkg.test(; coverage=true)'

.PHONY: init test coverage update
