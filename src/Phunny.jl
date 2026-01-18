"""
A semi-classical approach to solve for phonon frequencies and eigenmodes, built for integration with Sunny.jl.
Extension of Sunny.jl to include phonon analysis.
"""
module Phunny

using LinearAlgebra, SparseArrays, StaticArrays

#------------#
# Public API #
#------------#

export Model, Bond, build_model, neighbor_bonds_cutoff, neighbor_bonds_from_sunny
export assemble_force_constants!, assign_force_constants!, enforce_asr!

export dynamical_matrix, dynamical_gradient!, dynamical_hessian!
export phonons, group_velocities, msd_from_phonons, B_isotropic_from_phonons, U_from_phonons
export onephonon_dsf, onephonon_dsf_4d

export make_physical!, longitudinal_weights!
export ω_grid, mass_vector, q_cartesian, collapse, mass_lookup, bcoh_lookup, atomic_index

include("constants.jl")
include("helpers.jl")
include("calculations.jl")

#----------------------------#
# Validation & sanity checks #
#----------------------------#

export list_functions, validate_summary
export asr_residual, rigid_translation_residual, rigid_rotation_residual
export gamma_acoustic_energies, min_eigen_energy_along
       
include("devtests.jl")
include("systests.jl")


#-------------------
# Developer API
#-------------------

export AutoTargets, fit_per_pair!

include("ForceConstantSearch.jl")

end # module Phunny
