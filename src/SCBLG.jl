module SCBLG

    using Parameters
    using Quantica
    using ProgressMeter
    using LinearAlgebra
    using FixedPoint
    using Statistics
    using ArnoldiMethod
    using LinearMaps
    using SkewLinearAlgebra
    using Distributed
    using SIAMFANLEquations

    -#! from builderls.jl
    export build_meanfield_bilayergraphene
    export Yukawa, cutoff # finite-range potential 

    #! from selfconsistency.jl 
    # main tools for selfconsistency
    export build_meanfield_selfenergies # builds meanfield
    export calculate_fixedpoint # executes self-consistency
    export rasterscan_fixedpoint # parallelization for phase-digram
    export converged_hMF # builds hamiltonian from solution

    #! from analysis.jl
    # tools to analyse the results

    #! for testing from wherever
    export build_meanfield_selfenergies_onsite

    include("builders.jl")
    include("selfconsistency.jl")
    include("analysis.jl")

end
