The code foundations are written within the package Quantica.jl. 

Folder structure:

- /src : Julia source code

    - /src/builders : builds the bylayer graphene Hamiltonian
        Builds the lattice and tight-binding model and, from them, the self-consistent meanfield Hamiltonian w/ finite-range interactions.

    - /src/selfconsistency : calculates the self-consistent solutions and builds the meanfield Hamiltonian. 
        Builds the mean-field self-consistent iteration scheme by calculating the reduced density matrix ρ{i} of the previously built Hamiltonian using ρ{i-1} for the Hartree and Fock selfenergies Σ(ρ), starting from a hand-picked seed ρ{0}.

    - /src/analysis : tools to analyze the meanfield Hamiltonians.
        ...

    - /src/SCBLG.jl : module compiling the above

- /simulation : simulations to calculate the converged meanfield Hamiltonians.

    The simulation raw-data is large and must be generated using the simulations scripts prior to plotting the figures

- /analysis : calculation and plotting of results from the converged meanfield Hamiltonians.

    Even the figure's data is too large so it must also be generated using the simulated raw-data and the scripts at the top of each fig.jl file.