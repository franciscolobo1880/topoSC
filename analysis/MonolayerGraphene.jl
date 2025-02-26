#TODO

# from module
using Revise
using topoSC
# for specific file
using Quantica
using GLMakie

# symmetry points
const K = 2π/3 .* SA[+1, -1]
const M = 2π .* SA[1, 0] ./ 2

# visualize the hamiltonian
LatticeParams = (; ϵ=0.5)
TuningParams = (; μ=-3.0)
ph = build_monolayer_graphene(; LatticeParams...) # parametric
h = ph(; TuningParams...)
qplot(h; siteradius=0.05)

# visualize band structure along the symmetry path
LatticeParams = (; ϵ=0.0)
TuningParams = (; μ=-0.0)
plot_band_structure_KΓK(LatticeParams, TuningParams)

# visualize band structure near Dirac point
LatticeParams = (; ϵ=0.0)
TuningParams = (; μ=0.0)
plot_band_structure_K(LatticeParams, TuningParams)

#* ----- customizable plot functions -----

function plot_band_structure_KΓK(LatticeParams, TuningParams)
    ph = build_monolayer_graphene(; LatticeParams...) # parametric
    h = ph(; TuningParams...)
    b = bands(h, subdiv((0, 1, 3, 5, 6), 1001); mapping = (0, 1, 3, 5, 6) => (-M, -K, :Γ, +K, +M))
    
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1,1])
    qplot!(b; hide=(:nodes))

    ax.xlabel = L"ϕa_0"
    ax.xlabelsize = 30
    ax.xticklabelsize = 20
    ax.xticks = ([1, 3, 5], [L"K_-", L"Γ", L"K_+"])
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticklabelsize = 20

    save("figures/MonolayerGraphene/MonolayerGraphene_band_structure_KΓK.png", fig)
    return display(fig)
end

function plot_band_structure_K(LatticeParams, TuningParams)
    ph = build_monolayer_graphene(; LatticeParams...) # parametric
    h = ph(; TuningParams...)
    b = bands(h, subdiv((0, 1, 2), 1001); mapping = (0, 1, 2) => ( K-0.1*K, K, K-0.1*(M-K) ))
    
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1,1])
    qplot!(b; hide=(:nodes))

    ax.xlabel = L"ϕ a_0"
    ax.xlabelsize = 30
    ax.xticklabelsize = 20
    ax.xticks = ([1], [L"K_+"])
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticklabelsize = 20
    ylims!(-1.0, +1.0)

    save("figures/MonolayerGraphene/MonolayerGraphene_band_structure_K.png", fig)
    return display(fig)
end