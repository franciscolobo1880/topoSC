using Revise
using SCBLG
using GLMakie
using Quantica

import SCBLG.K as K
import SCBLG.δK as δK
const path = "/home/lobo/juliadev/SCBLG/analysis/localtesting/"

function calculating_lowenergy_bands(; Ez=0.0, μ=0.0, U=0.0)
    params = (; μ=μ, Ez=Ez, U=U)
    h = build_meanfield_bilayergraphene(nV=0.0)(; params...)
    b = bands(h, subdiv((0, 1, 2), 1001); mapping = (0, 1, 2) => (K-δK, K, K+δK))
    return b
end

# tweaking Ez with μ=0 and U=0 
begin
    Ezs = range(0.0, 0.01, 6)
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1,1], xticks=([0, 1, 2], [L"K-δK", L"K", L"K+δK"]), xlabel = L"ϕa_0", ylabel = L"ε \text{ (eV)}", xlabelsize=20, ylabelsize=20, titlesize=20, xticklabelsize=20, yticklabelsize=20)
    record(fig, path*"BLG_bands_versus_μ.gif", enumerate(Ezs); framerate=1) do (iEz, Ez)
        b = calculating_lowenergy_bands(; Ez)
        empty!(ax)   
        qplot!(b, hide=(:nodes), color=:black)
        ax.title = L"$\text{Low-energy bandstructure for }E_z=%$Ez$"
        #xlims!(K[1]+δK[1], K[1]+δK[2])
        ylims!(-0.02, +0.02)
    end
end