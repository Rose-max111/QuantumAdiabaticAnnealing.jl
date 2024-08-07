module QuantumAdiabaticAnnealing

using UnitDiskMapping
using Graphs
using KrylovKit
using CairoMakie
using Random
using CUDA
using GenericTensorNetworks
using SparseArrays
using DormandPrince
using LinearAlgebra
using Enzyme
using BitBasis

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
export track_equilibration_gausspulse_cpu!, SimulatedAnnealingHamiltonian, track_equilibration_gausspulse_gpu!

include("point.jl")
include("ca1d.jl")
include("quantum/energy_calculate.jl")
include("gadgets/gadgets.jl")
include("toy_model.jl")
include("cusa.jl")
include("autodiff.jl")
include("spinglass_adiabatic.jl")
include("visualize.jl")

end
