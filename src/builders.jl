#* ----- Common parameters -----

# pauli matrices: σ ≡ spin, τ ≡ nambu
const σ0 = τ0 = @SMatrix[1 0; 0 1]
const σx = τx = @SMatrix[0 1; 1 0]
const σy = τy = @SMatrix[0 -im; im 0]
const σz = τz = @SMatrix[1 0; 0 -1]
const ħ2me = 76.2 # meV h²/mₑ

#* ----- Kitaev chain -----
@with_kw struct LatticeParams_Kitaev @deftype Float64
    L  = Inf # chain length (nm)
    a0 = 1.0 #  lattice spacing (nm)
end;

build_Kitaev_chain(; kw...) = build_Kitaev_chain(LatticeParams_Kitaev(; kw...))
function build_Kitaev_chain(p::LatticeParams_Kitaev)
    @unpack L, a0 = p

    lat = LP.linear(; a0)

    model_normal = -@onsite((; μ=0.0) -> μ*τz) - @hopping((; t=0.0) -> t*τz)
    model_anomalous = @hopping((; Δ=0.0) -> Δ*1im*τy, region = (r, dr) -> dr[1] > 0)
    model_anomalous_dagger = @hopping((; Δ=0.0) -> -Δ*1im*τy, region = (r, dr) -> dr[1] < 0)
    model = model_normal + model_anomalous + model_anomalous_dagger

    h = lat |> hamiltonian(model, orbitals=2)

    if isfinite(L) h=supercell(h, region = r -> 0 <= r[1] <= L) end

    return h
end

#* ----- Oreg-Lutchyn nanowire -----
@with_kw struct LatticeParams_OregLutchyn @deftype Float64
    L  = Inf    # nanowire length (nm)
    a0 = 10.0   # lattice spacing (nm)
    m0 = 0.023  # effect mass of InAs
end;

build_OregLutchyn_wire(; kw...) = build_OregLutchyn_wire(LatticeParams_OregLutchyn(; kw...))
function build_OregLutchyn_wire(p::LatticeParams_OregLutchyn)
    @unpack L, a0, m0 = p
    t = ħ2me/(2m0*a0^2) # (meV)

    lat = LP.linear(; a0)

    model_K = @onsite((; μ=0.0) -> (2*t-μ)*kron(τz, σ0)) - hopping(t*kron(τz, σ0))
    model_Z = @onsite((; Vz=0.0) -> Vz*kron(τz, σz))
    model_α = @hopping((r, dr; α=0.0) -> α*(im*dr[1]/(2a0^2))*kron(τz, σy))
    model_Δ = @onsite((; Δ=0.0) -> Δ*kron(τy, σy))
    model = model_K + model_Z + model_α + model_Δ

    h = lat |> hamiltonian(model, orbitals=4)

    if isfinite(L) h=supercell(h, region = r -> 0 <= r[1] <= L) end

    return h
end

#* ----- Monolayer Graphene (MLG) -----
@with_kw struct LatticeParams_MLG @deftype Float64
    a0 = 0.246 # lattice spacing (nm)
    ϵ = 0.5 # onsite energy difference (eV)
    t = 3.16 # intralayer hopping strenght(eV)
end

build_monolayer_graphene(; kw...) = build_monolayer_graphene(LatticeParams_MLG(; kw...))
function build_monolayer_graphene(p::LatticeParams_MLG)
    @unpack a0, ϵ, t = p

    lat = LP.honeycomb(; a0, names = (:A, :B))

    model_A = @onsite((; μ=0.0) -> ϵ-μ, sublats=:A)
    model_B = @onsite((; μ=0.0) -> -ϵ-μ, sublats=:B)
    model_t = -hopping(t; range=a0/√3, sublats = (:B, :A) .=> (:A, :B))
    model = model_A + model_B + model_t
    
    h = hamiltonian(lat, model)

    return h
end

#* ----- Bilayer Graphene (BLG) -----
@with_kw struct LatticeParams_BLG @deftype Float64
    a0 = 0.246 # lattice spacing (nm)
    az = 1.36*a0 # layer separation length (nm)
    ϵD = 0.0 # D onsite energy (eV)
    ϵA = ϵD # A onsite energy (eV)
    ϵB = 0.022 # B onsite energy (eV)
    ϵC = ϵB # C onsite energy (eV)
    γ0 = 3.16 # intralayer hopping strenght C <-> D, A <-> B (eV)
    γ1 = 0.381 # interlayer hopping strenght B <-> C (eV)
    γ3 = 0.38 # interlayer hopping strenght A <-> C (eV)
    γ4 = 0.14 # interlayer hopping strenght A <-> C (eV)
end

build_bilayer_graphene(; kw...) = build_bilayer_graphene(LatticeParams_BLG(; kw...))
function build_bilayer_graphene(p::LatticeParams_BLG)
    @unpack a0, az, ϵD, ϵC, ϵB, ϵA, γ0, γ1, γ3, γ4 = p

    lat_top = LP.honeycomb(; a0, names = (:A, :B), dim=3) |> translate([0, 0, +az/2])
    lat_bot = LP.honeycomb(; a0, names = (:C, :D), dim=3) |> translate([0, a0/√3, -az/2])
    lat = combine(lat_bot, lat_top)

    model_top = @onsite((; μ=0.0) -> ϵB-μ, sublats=:B) + @onsite((; μ=0.0) -> ϵA-μ, sublats=:A) - hopping(γ0; range=a0/√3, sublats = (:B, :A) .=> (:A, :B))
    model_bot = @onsite((; μ=0.0) -> ϵD-μ, sublats=:D) + @onsite((; μ=0.0) -> ϵC-μ, sublats=:C) - hopping(γ0; range=a0/√3, sublats = (:D, :C) .=> (:C, :D))
    model_j1 = hopping(γ1; range=az, sublats = (:B, :C) .=> (:C, :B))
    model_j3 = -hopping(γ3; range=hypot(a0/√3, az), sublats = (:D, :A) .=> (:A, :D))
    model_j4 = hopping(γ4; range=hypot(a0/√3, az), sublats = (:D, :B) .=> (:B, :D)) + hopping(γ4; range=hypot(a0/√3, az), sublats = (:C, :A) .=> (:A, :C))
    model_E = @onsite((r; E=0.0) -> E*r[3])
    model = model_top + model_bot + model_j1 + model_j3 + model_j4 + model_E
    
    h = hamiltonian(lat, model)

    return h
end

