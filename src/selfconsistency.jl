#! main tools for self-consistency -----------------------------------------------------------------------------------------

# calculates the reduced densitry matrix ρ through g and then, from it, calculates the Hartree and Fock self-energies Σ(ρ)
function build_meanfield_selfenergies_onsite(g, V; params...)
    mf = meanfield(g; charge=kron(SCBLG.τz, SCBLG.σ0), nambu=true, namburotation=false, 
    onsite=params[:U])
    return mf
end

function build_meanfield_selfenergies(g, V; params...)
    mf = meanfield(g; charge=kron(SCBLG.τz, SCBLG.σ0), nambu=true, namburotation=false, 
    potential=r->V(r), onsite=params[:U], selector=(; range=params[:nV]*params[:a0]))
    return mf
end

# executes the iterative self-consistency routine with initial condition x0 from Anderson acceleration method
function calculate_fixedpoint(ΣMF::Quantica.MeanField; m=1, maxit=100, beta=0.8, atol=1e-7, ΔSC=2.0, ΔSM=0.3, params...)
    Σ0 = ΣMF(μ=0.0, kBT=0.0; params..., ΔSC, ΔSM) # acts a serializer translator 
    x0 = serialize(Float64, Σ0)

    function f!(x, x0, (ΣMF, Σ0, params))
        Σnew = Quantica.call!(ΣMF; params..., Σ=deserialize(Σ0, x0))
        copy!(x, serialize(Float64, Σnew))
        return x
    end

    Vstore = similar(x0, length(x0), 3m+3)
    sol = aasol(f!, x0, m, Vstore; pdata = (ΣMF, Σ0, params), maxit, beta, atol)
    sol.errcode == 0 || println("`aasol` unsuccessful at $params")
    return (; cΣ=sol.solution, i=length(sol.history), e=last(sol.history))
    # Σ≡converged Σ, i≡iterations and e≡error
end

# parallelization for phase-diagrams
function rasterscan_fixedpoint(fixedpoint, Vzs, μs)
    sols = @showprogress pmap(xy -> fixedpoint(xy[1], xy[2]), Iterators.product(Vzs, μs))
    cΣs = (x -> x.cΣ).(sols)
    is = (x -> x.i).(sols)
    es = (x -> x.e).(sols)
    return cΣs, is, es
end

# builds the meanfield hamiltonian from the converged solution Σ
# a given Σ0 translator works for all (Vz, μ)
function converged_hMF(ΣMF::Quantica.MeanField, h, sols, Vzs, μs; params...)
    Σ0 = ΣMF(μ=0.0, kBT=0.0; params..., ΔSC=2.0, ΔSM=0.4)
    hMF(iVz, iμ; params2...) = h(; params..., Vz=Vzs[iVz], µ=µs[iμ], Σ=deserialize(Σ0, sols[1][iVz, iμ]), params2...)
    return hMF
end