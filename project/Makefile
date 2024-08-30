JL = julia --project=..

device ?= 1
width ?= 12
depth ?= 30
gauss_width ?= 1.0
gradient ?= 1.5

calculate:
	$(JL) toymodel_genearalannealing.jl $(width) $(depth) $(gauss_width) $(gradient) $(device)

.PHONY: calculate
