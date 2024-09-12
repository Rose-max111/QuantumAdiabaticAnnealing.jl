JL = julia --project=..

device ?= 1
width ?= 14
depth ?= 15
gauss_width ?= 1.0
gradient ?= 1.24

calculate:
	$(JL) toymodel_genearalannealing.jl $(width) $(depth) $(gauss_width) $(gradient) $(device)

.PHONY: calculate
