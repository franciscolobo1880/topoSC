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
LatticeParams = (; ϵD=0.0, ϵA=0.0, ϵB=0.022, ϵC=0.022)
TuningParams = (; μ=0.0, E=0.0)
ph = build_bilayer_graphene(; LatticeParams...) # parametric
h = ph(; TuningParams...)
qplot(h; siteradius=0.05, hide=(:bravais))

# visualize band structure along the symmetry path
LatticeParams = (; ϵD=0.0, ϵA=0.0, ϵB=0.022, ϵC=0.022)
TuningParams = (; μ=0.0, E=0.0)
plot_band_structure_KΓK(LatticeParams, TuningParams)

# visualize band structure near Dirac point
LatticeParams = (; ϵD=0.0, ϵA=0.0, ϵB=0.022, ϵC=0.022)
TuningParams = (; μ=0.0, E=0.0)
plot_band_structure_K(LatticeParams, TuningParams)

# visualize trigonal warping near Dirac point
LatticeParams = (; ϵD=0.0, ϵA=0.0, ϵB=0.022, ϵC=0.022)
TuningParams = (; μ=0.0, E=0.0)
plot_trigonalwarping(LatticeParams, TuningParams)

#* ----- customizable plot functions -----

function plot_band_structure_KΓK(LatticeParams, TuningParams)
    ph = build_bilayer_graphene(; LatticeParams...) # parametric
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

    save("figures/BilayerGraphene/BilayerGraphene_band_structure_KΓK.png", fig)
    return display(fig)
end

function plot_band_structure_K(LatticeParams, TuningParams)
    ph = build_bilayer_graphene(; LatticeParams...) # parametric
    h = ph(; TuningParams...)
    b = bands(h, subdiv((0, 1, 2), 1001); mapping = (0, 1, 2) => ( K-0.1*K, K, K-0.1*(M-K) ))
    
    fig = Figure(size=(8band_structure_K00, 800))
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

    save("figures/BilayerGraphene/BilayerGraphene_band_structure_K.png", fig)
    return display(fig)
end

function plot_trigonalwarping(LatticeParams, TuningParams; nϕ=301, δK=0.05, ϵmin=0.0, ϵmax=0.005, Δϵ=0.0001)
    ph = build_bilayer_graphene(; LatticeParams...) # parametric
    h = ph(; TuningParams...)
    ϕ1s = range(K[1]-δK, K[1]+δK, length = nϕ)
    ϕ2s = range(K[2]-δK, K[2]+δK, length = nϕ)
    b = bands(h, ϕ1s, ϕ2s)

    ϵs = range(ϵmin, ϵmax; step=Δϵ)
    colors = cgrad(:gist_rainbow, length(ϵs), categorical=true);
    
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1,1])
    Colorbar(fig[1,2], colormap=colors, limits=(ϵmin, ϵmax))
    qplot!(b[(:, :, 0.0)], hide = (:nodes, :wireframe), color=:red)
    for (iϵ, ϵ) in enumerate(ϵs[begin+1:end])
        qplot!(b[(:, :, ϵ)], hide = (:nodes, :wireframe), color=colors[iϵ])
    end
    vl = vlines!(ax, [K[1]])
    hl = hlines!(ax, [K[2]])

    # plot axis customization
    ax.title = L"\text{Energy countours around Dirac point (}∼ \text{meV)}"
    ax.titlesize = 25
    ax.xlabel = L"ϕ_1 a_0"
    ax.xlabelsize = 30
    ax.xticklabelsize = 20
    ax.ylabel = L"ϕ_2 a_0"
    ax.ylabelsize = 30
    ax.yticklabelsize = 20
    vl.color = :black
    vl.linestyle = :dash
    hl.color = :black
    hl.linestyle = :dash

    save("figures/BilayerGraphene/BilayerGraphene_trigonalwarping.png", fig)
    return display(fig)
end