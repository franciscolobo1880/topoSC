# pauli matrices: σ ≡ spin, τ ≡ nambu, Π ≡ valley
const σ0 = Π0 = τ0 = I = @SMatrix[1 0; 0 1]
const σx = Πx = τx = @SMatrix[0 1; 1 0]
const σy = Πy = τy = @SMatrix[0 -im; im 0]
const σz = Πz = τz = @SMatrix[1 0; 0 -1]
# Dirac points in BZ basis
const K = SA[+2π/3, -2π/3]
const δK = SA[-0.02, +0.02]

# lattice construction parameters
@with_kw struct LatticeParams @deftype Float64
    a0 = 0.246 # lattice spacing (nm)
    az = 1.36*a0 # (nm)
    nV = 0 # neighbors max range
    ϵD = 0.0 # D onsite energy (eV)
    ϵA = ϵD # A onsite energy (eV)
    ϵB = 0.022 # B onsite energy (eV)
    ϵC = ϵB # C onsite energy (eV)
    γ0 = 3.16 # intralayer hopping strenght C <-> D, A <-> B (eV)
    γ1 = 0.381 # interlayer hopping strenght B <-> C (eV)
    γ3 = 0.38 # interlayer hopping strenght A <-> C (eV)
    γ4 = 0.14 # interlayer hopping strenght A <-> C (eV)
end

build_meanfield_bilayergraphene(; kw...) = build_meanfield_bilayergraphene(LatticeParams(; kw...))
function build_meanfield_bilayergraphene(p::LatticeParams)
    @unpack a0, az, nV, ϵD, ϵC, ϵB, ϵA, γ0, γ1, γ3, γ4 = p

    # building bilayer graphene lattice
    lat_top = LP.honeycomb(; a0, names = (:A, :B), dim=3) |> translate([0, 0, az])
    lat_bot = LP.honeycomb(; a0, names = (:C, :D), dim=3) |> translate([0, a0/√3, 0])
    lat = combine(lat_bot, lat_top) |> supercell(Int(nV+1), 2)

    # building tight-binding model (top + bottom + join): kinetics + zeeman
    model_t = @onsite((; μ=0.0) -> (ϵB-μ)*kron(τz, σ0), sublats=:B) + @onsite((; μ=0.0) -> (ϵA-μ)*kron(τz, σ0), sublats=:A) - hopping(γ0*kron(τz, σ0); range=a0/√3, sublats = (:B, :A) .=> (:A, :B))
    model_b = @onsite((; μ=0.0) -> (ϵD-μ)*kron(τz, σ0), sublats=:D) + @onsite((; μ=0.0) -> (ϵC-μ)*kron(τz, σ0), sublats=:C) - hopping(γ0*kron(τz, σ0); range=a0/√3, sublats = (:D, :C) .=> (:C, :D))
    model_j1 = hopping(γ1*kron(τz, σ0); range=az, sublats = (:B, :C) .=> (:C, :B))
    model_j3 = -hopping(γ3*kron(τz, σ0); range=hypot(a0/√3, az), sublats = (:D, :A) .=> (:A, :D))
    model_j4 = hopping(γ4*kron(τz, σ0); range=hypot(a0/√3, az), sublats = (:D, :B) .=> (:B, :D)) + hopping(γ4*kron(τz, σ0); range=hypot(a0/√3, az), sublats = (:C, :A) .=> (:A, :C))
    model_Ez = @onsite((r; Ez=0.0) -> Ez*r[3]*kron(τz, σz))
    model_Δ = @onsite((; Δ=0.0) -> Δ*kron(τy, σy)) # [essencial that ΔSM = 0, just for structure]
    model_0 = model_t + model_b + model_j1 + model_j3 + model_j4 + model_Ez + model_Δ

    # hartree and fock meanfield self-energies w/ finite-range interactions
    ΣHartree = @onsite((i; Σ=zerofield) --> Σ[i]) 
    ΣFock = @hopping((i, j; Σ=zerofield) --> Σ[i, j]; range=nV*a0) 
    model = model_0 + ΣHartree + ΣFock

    # building (and clipping) hamiltonian 
    h = hamiltonian(lat, model; orbitals=4)

    return h
end

# yukawa potential
Yukawa(r; params, A, p::LatticeParams=LatticeParams()) = params.U*A*exp(-norm(r/p.a0)/params.λ)/norm(r/p.a0)
Yukawa(; params) = r -> Yukawa(r; params, A=1/2)  # Returns a function of r
# cutoff function [nV = 0 is onsite]
function cutoff(V; params, pU=1, p::LatticeParams=LatticeParams())
    @unpack a0 = p
    if iszero(params.λ) || isnothing(V)
        return 0
    else
        cut = abs(pU/100 * params.U / log10(params.λ))
        v = [n-1 for n in 2:1e3 if abs(V(n*a0)) < cut]
        return Int(first(v))
    end
end

