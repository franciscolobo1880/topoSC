#TODO topological invariant (bandstructure title changing accordingly)

# from module
using Revise
using topoSC
# for specific file
using Quantica
using GLMakie

# visualize the hamiltonian
params = (; μ=-3.0, t=1.0, Δ=1.0)
ph = build_Kitaev_chain(L=Inf)
h = ph(; params...)
qplot(h)

# visualize band structure
ks = subdiv(-π, +π, 1001)
plot_band_band_structure(h, ks)

# visualize band spectrum versus μ
μs = range(-3.0, 0.0; step=0.01)
plot_bands_gif(h, μs)

# band spectrum versus μ
μs = range(-3.0, 0.0; step=0.01)
ph = build_Kitaev_chain(L=5.0) 
b = bands(ph, μs; mapping = μ -> ftuple(;params..., Δ=1.0, t=1.0, μ), solver=ES.ShiftInvert(ES.ArnoldiMethod(nev = 1), 0.0001))
begin
    fig = Figure(size=(2400, 800))
    ax = Axis(fig[1, 1],)
    qplot!(b, color=:black, hide=(:nodes))

    ax.xlabel = L"μ\text{ (meV)}"
    ax.xlabelsize = 30
    ax.ylabel = L"ϵ\text{ (meV)}"
    ax.ylabelsize = 30
    xlims!(μs[begin], μs[end])
    ylims!(-2, 2)
    display(fig)
end

#* ----- prettier plots -----

function plot_band_structure(h, ks)
    b = bands(h, ks)

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1])
    qplot!(b, hide=(:bands), color=:red)
    ax.title = L"$μ=%$(params.μ),\ t=%$(params.t),\ Δ=%$(params.Δ)$"
    ax.titlesize = 30
    ax.xlabel = L"k a_0"
    ax.xlabelsize = 30
    ax.xticks=([-π, -π/2, 0, π/2, π], [L"-π", L"-π/2", L"Γ", L"π/2", L"π"])
    ax.xticklabelsize = 20
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticks = -6:6 
    ylims!(-6.0, +6.0)

    save(folder*"figures/Kitaev_bands.png", fig)
    return display(fig)
end

function plot_bands_gif(h, μs)
    ks = subdiv(-π, π, 1000)
    fig = Figure(size=(600,500))
    ax = Axis(fig[1,1], ylabel=L"ϵ\text{ (meV)}", xlabel=L"k \text{ (nm)}^{-1}", xlabelsize=20, ylabelsize=20, xticks=([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]), yticks=(-60:10:60, string.(-60:10:60)))
    record(fig, "figures/Kitaev.gif", eachindex(μs); framerate=4) do iμ
        b = bands(ph(; params..., μ=μs[iμ]), ks)
        empty!(ax)   
        qplot!(b, color=:black, hide=(:bands))
    end
    return display(fig)
end


