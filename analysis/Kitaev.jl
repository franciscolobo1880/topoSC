#TODO topological invariant (bandstructure title changing accordingly)
#TODO Cesar's Majorana oscillations

# from module
using Revise
using topoSC
# for specific file
using Quantica
using GLMakie

# visualize the hamiltonian
LatticeParams = (; L=4, a0=1.0)
TuningParams = (; μ=-3.0, t=1.0, Δ=1.0)
ph = build_Kitaev_chain(; LatticeParams...) # paramatric
h = ph(; TuningParams...)
qplot(h)

# visualize band structure
LatticeParams = (; L=Inf, a0=1.0)
TuningParams = (; μ=-3.0, t=1.0, Δ=1.0)
ks = subdiv(-π, +π, 1001)
plot_band_structure(LatticeParams, TuningParams, ks)

# visualize band structure evolution versus μ
LatticeParams = (; L=Inf, a0=1.0)
TuningParams = (; μ=-3.0, t=1.0, Δ=1.0)
ks = subdiv(-π, +π, 1001)
μs = range(-3.0, 0.0; step=0.25)
plot_band_structure_gif(LatticeParams, TuningParams, ks, μs)

# visualize band spectrum versus μ
LatticeParams = (; L=50.0, a0=1.0)
TuningParams = (; μ=-3.0, t=1.0, Δ=1.0)
μs = range(-3.0, 0.0; step=0.01)
plot_band_spectrum(LatticeParams, TuningParams, μs)

#* ----- customizable plot functions -----

function plot_band_structure(LatticeParams, TuningParams, ks)
    @assert LatticeParams.L == Inf "Chain length needs to be infinite!"
    ph = build_Kitaev_chain(; LatticeParams...)
    h = ph(; TuningParams...)
    b = bands(h, ks)

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1])
    qplot!(b, hide=(:bands), color=:red)
    ax.title = L"$μ=%$(TuningParams.μ),\ t=%$(TuningParams.t),\ Δ=%$(TuningParams.Δ)$"
    ax.titlesize = 30
    ax.xlabel = L"k a_0 (nm^{-1})"
    ax.xlabelsize = 30
    ax.xticks=([-π, -π/2, 0, π/2, π], [L"-π", L"-π/2", L"Γ", L"π/2", L"π"])
    ax.xticklabelsize = 20
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticks = -6:6 
    ylims!(-6.0, +6.0)

    save("figures/Kitaev/Kitaev_band_structure.png", fig)
    return display(fig)
end

function plot_band_structure_gif(LatticeParams, TuningParams, ks, μs)
    @assert LatticeParams.L == Inf "Chain length needs to be infinite!"
    ph = build_Kitaev_chain(; LatticeParams...)

    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    ax.titlesize = 30
    ax.xlabel = L"k a_0 (nm^{-1})"
    ax.xlabelsize = 30
    ax.xticks=([-π, -π/2, 0, π/2, π], [L"-π", L"-π/2", L"Γ", L"π/2", L"π"])
    ax.xticklabelsize = 20
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticks = -6.0:6.0 
    ylims!(-6.0, +6.0)

    record(fig, "figures/Kitaev/Kitaev_band_structure.gif", enumerate(μs); framerate=4) do (iμ, μ)
        b = bands(ph(; TuningParams..., μ=μ), ks)
        empty!(ax)   
        qplot!(b, color=:black, hide=(:bands))
        ax.title = L"$t=%$(TuningParams.t),\ Δ=%$(TuningParams.Δ), \ μ=%$μ$"
    end

    return display(fig)
end

function plot_band_spectrum(LatticeParams, TuningParams, μs)
    @assert LatticeParams.L != Inf "Chain length needs to be finite!"
    ph = build_Kitaev_chain(; LatticeParams...)
    slvr = ES.ShiftInvert(ES.ArnoldiMethod(nev = Int(LatticeParams.L)), 0.0001)
    b = bands(ph, μs; mapping = μ -> ftuple(;TuningParams..., μ), solver=slvr)

    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1],)
    qplot!(b, color=:black, hide=(:nodes))

    ax.xlabel = L"μ\text{ (meV)}"
    ax.xlabelsize = 30
    ax.ylabel = L"ϵ\text{ (meV)}"
    ax.ylabelsize = 30
    xlims!(μs[begin], μs[end])
    ylims!(-2, 2)

    save("figures/Kitaev/Kitaev_band_spectrum.png", fig)
    return display(fig)
end