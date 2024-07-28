JL = julia --project=..

device ?= 1
width ?= 10
depth ?= 10

calculate:
	$(JL) toymodel_hamiltonian_cuda.jl $(width) $(depth) $(gradient) $(device)

.PHONY: calculate
