JL = julia --project=..

device ?= 1

calculate:
	$(JL) toymodel_hamiltonian_cuda.jl $(width) $(depth) $(gradient) $(device)

.PHONY: calculate
