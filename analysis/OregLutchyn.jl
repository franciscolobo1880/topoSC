#TODO

# from module
using Revise
using topoSC
# for specific file
using Quantica
using GLMakie

# visualize the hamiltonian
LatticeParams = (; L=40.0, a0=10.0, m0=0.023)
TuningParams = (; α=40.0, Δ=0.3, μ=0.0, Vz=0.0)
ph = build_OregLutchyn_wire(; LatticeParams...) # parametric
h = ph(; TuningParams...)
qplot(h)

# visualize band structure
LatticeParams = (; L=Inf, a0=10.0, m0=0.023)
TuningParams = (; α=40.0, Δ=0.3, μ=0.0, Vz=0.0)
ks = subdiv(-π, +π, 1001)
plot_band_structure(LatticeParams, TuningParams, ks)

# visualize band structure evolution versus Zeeman
LatticeParams = (; L=Inf, a0=10.0, m0=0.023)
TuningParams = (; α=40.0, Δ=0.3, μ=0.0, Vz=0.0)
ks = subdiv(-π, +π, 1001)
Vzs = range(0.0, 1.0; step=0.1)
plot_band_structure_gif(LatticeParams, TuningParams, ks, Vzs)

# visualize band spectrum versus Zeeman
LatticeParams = (; L=2000.0, a0=10.0, m0=0.023)
TuningParams = (; α=40.0, Δ=0.3, μ=0.0, Vz=0.0)
Vzs = range(0.0, 1.0; step=0.01)
plot_band_spectrum(LatticeParams, TuningParams, Vzs)

# visualize phase diagram
LatticeParams = (; L=Inf, a0=10.0, m0=0.023)
TuningParams = (; α=40.0, Δ=0.3, μ=0.0, Vz=0.0)
Vzs = range(0.0, 1.0; step=0.01)
μs = range(-1.0, 1.0; step=0.01)
plot_phase_diagram(LatticeParams, TuningParams, Vzs, μs)

#* ----- customizable plot functions -----

function plot_band_structure(LatticeParams, TuningParams, ks)
    @assert LatticeParams.L == Inf "Nanowire length needs to be infinite!"
    ph = build_OregLutchyn_wire(; LatticeParams...)
    h = ph(; TuningParams...)
    b = bands(h, ks)

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1])
    qplot!(b, hide=(:nodes), color=:red)
    ax.title = L"$μ=%$(TuningParams.μ),\ α=%$(TuningParams.α),\ Δ=%$(TuningParams.Δ),\ V_z=%$(TuningParams.Vz)$"
    ax.titlesize = 30
    ax.xlabel = L"k a_0 (nm^{-1})"
    ax.xlabelsize = 30
    ax.xticks=([-π, -π/2, 0, π/2, π], [L"-π", L"-π/2", L"Γ", L"π/2", L"π"])
    ax.xticklabelsize = 20
    ax.ylabel = L"ε \text{ (meV)}"
    ax.ylabelsize = 30
    ax.yticks = -10.0:10.0 
    xlims!(-0.5, +0.5)
    ylims!(-2.0, +2.0)

    save("figures/OregLutchyn/OregLutchyn_band_structure.png", fig)
    return display(fig)
end

function plot_band_structure_gif(LatticeParams, TuningParams, ks, Vzs)
    @assert LatticeParams.L == Inf "Nanowire length needs to be infinite!"
    ph = build_OregLutchyn_wire(; LatticeParams...)

    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    ax.titlesize = 30
    ax.xlabel = L"k a_0 (nm^{-1})"
    ax.xlabelsize = 30
    ax.xticks=([-π, -π/2, 0, π/2, π], [L"-π", L"-π/2", L"Γ", L"π/2", L"π"])
    ax.xticklabelsize = 20
    ax.ylabel = L"ε \text{ (eV)}"
    ax.ylabelsize = 30
    ax.yticks = -10.0:10.0 
    xlims!(-0.5, +0.5)
    ylims!(-2.0, +2.0)

    record(fig, "figures/OregLutchyn/OregLutchyn_band_structure.gif", enumerate(Vzs); framerate=4) do (iVz, Vz)
        b = bands(ph(; TuningParams..., Vz=Vz), ks)
        empty!(ax)   
        qplot!(b, color=:black, hide=(:nodes))
        ax.title = L"$μ=%$(TuningParams.μ),\ α=%$(TuningParams.α),\ Δ=%$(TuningParams.Δ),\ V_z=%$(Vz)$"
    end

    return display(fig)
end

function plot_band_spectrum(LatticeParams, TuningParams, Vzs)
    @assert LatticeParams.L != Inf "Nanowire length needs to be finite!"
    ph = build_OregLutchyn_wire(; LatticeParams...)
    slvr = ES.ShiftInvert(ES.ArnoldiMethod(nev = 50), 0.0001)
    b = bands(ph, Vzs; mapping = Vz -> ftuple(;TuningParams..., Vz), solver=slvr)

    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1],)
    qplot!(b, color=:black, hide=(:nodes))

    ax.title = L"$μ=%$(TuningParams.μ),\ α=%$(TuningParams.α),\ Δ=%$(TuningParams.Δ)$"
    ax.titlesize = 30
    ax.xlabel = L"μ\text{ (meV)}"
    ax.xlabelsize = 30
    ax.ylabel = L"ϵ\text{ (meV)}"
    ax.ylabelsize = 30
    xlims!(Vzs[begin], Vzs[end])
    ylims!(-0.5, 0.5)

    save("figures/OregLutchyn/OregLutchyn_band_spectrum.png", fig)
    return display(fig)
end


function plot_phase_diagram(LatticeParams, TuningParams, Vzs, μs)
    @assert LatticeParams.L == Inf "Nanowire length needs to be infinite!"
    ph = build_OregLutchyn_wire(; LatticeParams...)
    Ωs = [Quantica.gap(ph(;TuningParams..., Vz=Vz, μ=μ); nev=2) for Vz in Vzs, μ in μs]
    Pfaffian(μ, Δ) = sqrt(μ^2+Δ^2)

    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, Vzs, μs, Ωs)    
    bar = Colorbar(fig[1, 2], hm)    
    Pfa = lines!(Pfaffian.(μs, TuningParams.Δ), μs) 

    ax.title = L"$α=%$(TuningParams.α),\ Δ=%$(TuningParams.Δ)$"
    ax.titlesize = 30
    ax.xlabel = L"V_\text{Z}\text{ (meV)}"
    ax.xlabelsize = 30
    ax.ylabel = L"μ\ \text{(meV)}"
    ax.ylabelsize = 30
    hm.colormap = cgrad(:inferno)
    hm.colorrange = (0,0.401)
    hm.highclip = cgrad(:inferno)[end]    
    bar.label = L"Ω\text{ (meV)}"
    bar.labelsize = 30
    bar.ticksvisible = false
    bar.ticklabelsvisible = false
    Pfa.color = :white
    Pfa.linestyle = :dash
    xlims!(Vzs[begin], Vzs[end])
    ylims!(μs[begin], μs[end])

    save("figures/OregLutchyn/OregLutchyn_phase_diagram.png", fig)
    return display(fig)
end