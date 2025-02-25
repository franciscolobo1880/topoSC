using Revise
using SCBLG
using GLMakie
using Quantica

params = (; μ=0.0, Ez=0.0, U=0.0)
h = build_meanfield_bilayergraphene(nV=0)(; params...)
qplot(h; siteradius=0.05, hide=(:bravais))

h2D = h |> supercell(4)
qplot(h2D; siteradius=0.05, hide=(:bravais))

h1D = h |> stitch((:, 0.1))
qplot(h1D; siteradius=0.05, hide=(:bravais))

h1D´ = @stitch(h, SA[1], ϕ)




