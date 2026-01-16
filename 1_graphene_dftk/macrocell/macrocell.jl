using DFTK, LinearAlgebra
using ASEconvert
include("../../lib/band_diagram.jl")
px = println

dist(a,b) = 2*abs(a-b)/(abs(a) + abs(b))

function scf_graphene()
    ecut = 15
    a = 4.66 # length of the vectors of the fundamental cell, in Bohr
    L = 20  # height of the simulation box
    kgrid = [6, 6, 1]
    # Loads model
    file_psp = "hgh/pbe/c-q4"
    psp = load_psp(file_psp)
    # Choose carbon atoms
    C = ElementPsp(:C, psp=psp)
    # Puts carbon atoms at the right positions for graphene
    C1 = [1/3,-1/3,0.0]; C2 = -C1

    a1 = a*[1/2,-sqrt(3)/2, 0]
    a2 = a*[1/2, sqrt(3)/2, 0]
    a3 = L*[0  , 0        , 1]
    lattice = [a1 a2 a3]

    atoms = [C,C]
    positions = [C1,C2]

    repeat = 2

    # Supercell
    unit_cell = periodic_system(lattice, atoms, positions)

    function solve_n(n)
        supercell_ase = convert_ase(unit_cell) * pytuple((n, 1, 1))
        supercell = pyconvert(AbstractSystem, supercell_ase)
        supercell = attach_psp(supercell; C=file_psp)

        # Builds model
        # model = model_PBE(lattice, atoms, positions; temperature=1e-3, smearing=Smearing.Gaussian())
        model = model_PBE(supercell; temperature=1e-3, smearing=Smearing.Gaussian())
        # Builds basis
        basis = DFTK.PlaneWaveBasis(model; Ecut=ecut, kgrid)
        # px("Number of spin components ",model.n_spin_components)
        # Loads sizes
        (N,N,Nz) = basis.fft_size
        px("For ",n," fftsize is ",N," ",Nz)
        scfres = self_consistent_field(basis)
        px("Energies")
        Ens = scfres.energies.values
        # px(Ens)
        E = sum(Ens)/n
        E
    end
    Es = []
    Nmax = 2
    for n=1:Nmax
        E = solve_n(n)
        push!(Es,E)
        px("n ",n," E=",E," dist ",dist(E,Es[1]))
    end
    dE = [dist(Es[i],Es[1]) for i=1:Nmax]
    px(dE)
end

scf_graphene()
