using Revise
using SCBLG
using GLMakie
using Quantica

import SCBLG.K as K
import SCBLG.δK as δK
path = "/home/lobo/juliadev/SCBLG/analysis/localtesting/"

#! defining parameters
params = (; U=0.1)
μs = range(0.1, 0.1, 1)
Ezs = range(0.1, 0.1, 1)

#! self-consistency iterations
h = build_meanfield_bilayergraphene(nV=0.0)
ΣMF = build_meanfield_selfenergies_onsite(greenfunction(h, GS.Schur()), nothing; params...)
fp(Ez, μ) = calculate_fixedpoint(ΣMF; params..., Ez, μ)
sols = rasterscan_fixedpoint(fp, Ezs, μs)
sols[2] # checking number of iterations 

#! construction meanfield hamiltonian from solution
Σ0 = ΣMF(μ=0.0, kBT=0.0; params...)
hMF(iEz, iμ; params...) = h(; params..., Ez=Ezs[iEz], µ=µs[iμ], Σ=deserialize(Σ0, sols[1][iEz, iμ]), params...)

#! analysis in a specific parameter combinations
hMFℓ = hMF(1, 1)

# full bandstructure
qplot(bands(h), hide=(:nodes, :wireframe))

# bands along symmetry path
b = bands(h, subdiv((0, 2, 3), 1001); mapping = (0, 2, 3) => (:Γ, :K, :M))
qplot(b, axis = (; xticks = ([0, 2, 3], [L"Γ", L"K_+", L"M"]), xlabel = L"ka_0", ylabel = L"ε \text{ (eV)}"))

# low-energy bands along symmetry path
b = bands(h, subdiv((0, 1, 2), 1001); mapping = (0, 1, 2) => (K-δK, K, K+δK))
qplot(b, axis = (; xticks = ([0, 1, 2], [L"K_+ - δK", L"K_+", L"K_+ + δK"]), xlabel = L"ka_0", ylabel = L"ε \text{ (eV)}"))

# low-energy bands tringular warping
ϕ1s = range(K[1]+δK[1], K[1]+δK[2], length = 301);
ϕ2s = range(K[2]+δK[1], K[2]+δK[2], length = 301);
b = bands(h, ϕ1s, ϕ2s)
qplot(b, hide=(:nodes, :wireframe))
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