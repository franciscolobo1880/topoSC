#!/usr/bin/env -S julia --project

## SBATCH --nodes=1
#SBATCH --partition=most
#SBATCH --ntasks=192
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output="intrinsic_infinite.out"

# Launch workers
using Distributed
using ProgressMeter
const maxprocs = 96
addprocs(max(0, maxprocs + 1 - nworkers()))

# Load packages
import Serialization
@everywhere begin
    using Quantica
    using HartreeFockWires
    using SIAMFANLEquations
    Quantica.BLAS.set_num_threads(1)
end

# ----- Self-consistency simulations -----
@everywhere params = (; U=-1.0, λ=0.0) 

@everywhere μs = range(-1.0, -1.0, 101) # meV (ħ²/mₑ)
@everywhere Ezs = range(0.0, 1.0, 101) # meV (ħ²/mₑ)

# Yukawa potential
@everywhere V(r) = Yukawa(r; params)
@everywhere params = merge(params, (; nV=cutoff(V; params)))

@everywhere begin
    Quantica.BLAS.set_num_threads(1)

    h = build_meanfield_bilayergraphene(U=params.U, nV=params.nV)
    ΣMF = build_meanfield_selfenergies(greenfunction(h, GS.Schur()), V; params...)
    fp(Ez, μ) = calculate_fixedpoint(ΣMF; params..., Ez, μ)
end

begin
    @time sols = rasterscan_fixedpoint(fp, Ezs, μs)
    path = "analysis/phasediagram/rawdata/BLG_$(params).dat" 
    Serialization.serialize(path, (params, Ezs, μs, sols))
end

# -- cleanup --
println("Finished!")
rmprocs(workers())