begin
    using Revise
    using SCBLG
    import Serialization
    using Quantica
    using ProgressMeter
    using CairoMakie
    CairoMakie.activate!(type = "png")
end;

#! loading raw data and calculating gap and spectrum----------------------------------------------------------------
# calculating phase diagram
begin
    path = "analysis/phasediagram/rawdata/BLG_(U = -1.0, nV = 2).dat"
    (params, Ezs, μs, sols) = Serialization.deserialize("$path")
    h = build_meanfield_bilayergraphene(U=params.U, nV=params.nV)
    ΣMF = build_meanfield_selfenergies(greenfunction(h, GS.Schur()), nothing; params...)
    hMF = converged_hMF(ΣMF, h, sols, Ezs, μs; params...)
    Ωs = [Quantica.gap(hMF(;params..., Ez=Ez, μ=μ); nev=2) for Ez in Ezs, μ in μs]
    data = (; params, Ezs, μs, hMF, Ωs)
    return data
end;

# picking by hand iμℓs from the phase diagram
begin
    fig = Figure(size=(900, 700))
    ax = Axis(fig[1, 1], xlabel=L"E_z \text{ (meV)}", ylabel = L"μ\text{ (meV)}")
    hm = heatmap!(ax, Ezs, μs, data.Ωs, colormap=cgrad(:inferno, scale=:exp), colorrange=(0.0, 0.4), highclip=cgrad(:inferno, scale=:exp)[end]) 
    Colorbar(fig[1, 2], hm, label=L"Ω\text{ (meV)}", labelsize=22.5)
    hlines!(ax, μs[30], linestyle=Linestyle([0, 15, 20]), color=:red)
    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    display(fig)
    save(folder*"fig1.pdf", fig)
end

# calculing the band spectrum at those iμℓs
iμℓs = [10, 20, 30]
begin
    BSs = Vector{Any}(undef, length(iμℓs))
    solver(Ei) = ES.ShiftInvert(ES.ArnoldiMethod(nev = Ei), 0)
    for (i, iμ) in enumerate(iμℓs)
        BSs[i] = bands(hMF, Ezs; mapping=Ez->ftuple(; params..., Ez, µ=µs[iμ], ΣMF), solver=solver(4))
    end

    data = merge(data, (; iμℓs=iμℓs))
    data = merge(data, (; BSs=BSs))
    path = "analysis/phasediagram/BLG_(U = -1.0, nV = 2).dat"
    Serialization.serialize(path, data)
end

#! loading if already calculated and saved ----------------------------------------------------------
folder = "analysis/phasediagram/BLG_(U = -1.0, nV = 2).dat"
data = Serialization.deserialize(folder*"fig2.dat")

#! plotting of the complete figure w/ details ------------------------------------------------------
begin
    fig = Figure(size=(900, 700))

    # plotting phase diagram versus μ
    ax = Axis(fig[1, 1], xlabel=L"E_z \text{ (meV)}", ylabel = L"μ\text{ (meV)}", xlabelsize=22.5, ylabelsize=22.5)
    hm = heatmap!(ax, Ezs, μs, data.Ω, colormap=cgrad(:inferno, scale=:exp), colorrange=(0.0, 0.4), highclip=cgrad(:inferno, scale=:exp)[end]) 
    Colorbar(fig[1, 2], hm, label=L"Ω\text{ (meV)}", labelsize=22.5)
    hlines!(ax, μs[data.iμℓs[1]], linestyle=Linestyle([0, 15, 20]), color=:red)
    #hlines!(ax, μs[data.iμℓs[2]], linestyle=Linestyle([0, 15, 20]), color=:limegreen)
    #hlines!(ax, μs[data.iμℓs[3]], linestyle=Linestyle([0, 15, 20]), color=:blue)
    text!(fig[1, 1], L"\text{s-wave metalic}", color=:black, fontsize=25, position=(0.1, -8.5))
    text!(fig[1, 1], L"\text{s-wave insulator}", color=:white, fontsize=25, position=(0.55, -9.0))
    text!(fig[1, 1], L"\text{p-wave metalic}", color=:black, fontsize=25, position=(0.1, -8.5))
    text!(fig[1, 1], L"\text{p-wave insulator}", color=:white, fontsize=25, position=(0.55, -9.0))

    # plotting band spectrums
    ax = Axis(fig[2, 1], xlabel=L"E_z \text{ (meV)}", ylabel = L"ϵ\text{ (meV)}", xlabelsize=22.5, ylabelsize=22.5, xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible=true, yminorticksvisible=true, yminorgridvisible=true)
    qplot!(data.BSs[1], color=:red, hide=:nodes)
    qplot!(data.BSs[2], color=:limegreen, hide=:nodes)
    qplot!(data.BSs[3], color=:blue, hide=:nodes)
    #vlines!(ax, Vzs[1], linestyle=:dot, color=:black) 

    rowsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    display(fig)
    save(folder*"fig1.pdf", fig)
end

