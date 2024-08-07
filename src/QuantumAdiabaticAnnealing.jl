module QuantumAdiabaticAnnealing

using BitBasis
using CairoMakie
using CUDA
using DormandPrince
using Enzyme
using GenericTensorNetworks
using Graphs
using KrylovKit
using LinearAlgebra
using Random
using SparseArrays

export Point
export generate_random_lattice
export mis_gadget, CellularAutomata1D, topology, logic_gate
export show_transversal_graph_weight
export track_equilibration!, SimulatedAnnealingMIS

export toy_model_state_energy, toy_model_transition_matrix
export TransitionRule, calculate_energy, local_energy, Metropolis, HeatBath

export sp_check_vaild, spinglass_adiabatic_dp8, sp_ground_state_sa, sp_ground_state, spinglass_mapping, sp_energy
export spinglassmodel, spinglass_mapping, instantaneous_field, spinglass_random_mapping

export spinglass_hamiltonian, instantaneous_field_autodiff, runge_kutta_integrate!

include("point.jl")
include("ca1d.jl")
include("gadgets/gadgets.jl")
include("toy_model.jl")
include("cusa.jl")
include("autodiff.jl")
include("spinglass_adiabatic.jl")
include("visualize.jl")

end
