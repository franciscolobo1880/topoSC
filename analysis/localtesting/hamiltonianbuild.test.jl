using Revise
using SCBLG
using GLMakie
using Quantica

import SCBLG.K as K
import SCBLG.δK as δK
const path = "/home/lobo/juliadev/SCBLG/analysis/localtesting/"

# hamiltonian
params = (; μ=0.0, Ez=0.0, U=0.0)
h = build_meanfield_bilayergraphene(nV=0)(; params...)
qplot(h; siteradius=0.05, hide=(:bravais))

qplot(bands(h), hide=(:nodes, :wireframe))

# bands along symmetry path
b = bands(h, subdiv((0, 2, 3), 1001); mapping = (0, 2, 3) => (:Γ, :K, :M))
qplot(b, axis = (; xticks = ([0, 2, 3], [L"Γ", L"K_+", L"M"]), xlabel = L"ka_0", ylabel = L"ε \text{ (eV)}"))

# low-enegy bands along symmetry path
b = bands(h, subdiv((0, 1, 2), 1001); mapping = (0, 1, 2) => (K-δK, K, K+δK))
qplot(b, axis = (; xticks = ([0, 1, 2], [L"K_+ - δK", L"K_+", L"K_+ + δK"]), xlabel = L"ka_0", ylabel = L"ε \text{ (eV)}"))

# low-energy bands tringular warping
ϕ1s = range(K[1]+δK[1], K[1]+δK[2], length = 301);
ϕ2s = range(K[2]+δK[1], K[2]+δK[2], length = 301);
b = bands(h, ϕ1s, ϕ2s)
#qplot(b, hide=(:nodes, :wireframe))
begin
    ϵmax=0.005; Δϵ=0.0001
    ϵs = range(0.0, ϵmax; step=Δϵ)
    colors = cgrad(:gist_rainbow, length(ϵs), categorical=true);
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1,1], xticks = ([K[1]+δK[1], K[1], K[1]+δK[2]], [L"K_{ϕ_1}-δK", L"K_{ϕ_1}", L"K_{ϕ_1}+δK"]), yticks=([K[2]+δK[1], K[2], K[2]+δK[2]], [L"K_{ϕ_2}-δK", L"K_{ϕ_2}", L"K_{ϕ_2}+δK"]), xlabel=L"ϕ_1 a_0", ylabel=L"ϕ_2 a_0", title=L"\text{Energy cuts in the order of meV}", xlabelsize=20, ylabelsize=20, titlesize=20, xticklabelsize=20, yticklabelsize=20)
    qplot!(b[(:, :, 0.0)], hide = (:nodes, :wireframe), color=:red)
    for (iϵ, ϵ) in enumerate(ϵs[begin+1:end])
        qplot!(b[(:, :, ϵ)], hide = (:nodes, :wireframe), color=colors[iϵ])
    end
    δKℓ = 0.002
    xlims!(K[1]+δK[1]-δKℓ, K[1]+δK[2]+δKℓ)
    ylims!(K[2]+δK[1]-δKℓ, K[2]+δK[2]+δKℓ)
    save(path*"BLG_triangularwarping.png", fig)
    display(fig)
end