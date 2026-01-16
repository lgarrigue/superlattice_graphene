# include("graphene.jl")
using JLD
using JLD2

# PATH_EXPORTS_MONOLAYER_root = "../1_graphene_dftk/apply_graphene_outputs/"
PATH_EXPORTS_MONOLAYER_root = "../apply_graphene_outputs/"
PATH_EXPORTS_MONOLAYER_root = "/home/louis/heavy_stuff/graphs_code/improved_graphene/apply_graphene_outputs/" # use a directory which is not synchronized with a cloud !
PATH_EXPORTS_MONOLAYER = PATH_EXPORTS_MONOLAYER_root*"exported_functions/"

function exports_v_u1_u2(p)
    create_dir(PATH_EXPORTS_MONOLAYER_root)
    create_dir(PATH_EXPORTS_MONOLAYER)
    filename_root = string(PATH_EXPORTS_MONOLAYER ,"N",p.basis.N,"_Nz",p.basis.Nz)
    filename = filename_root*".jld"
    save(filename,"N",p.basis.N,"Nz",p.basis.Nz,"a",p.basis.a,"L",p.basis.L,"v_dir",p.v_dir,"Es_K",p.Es_K,"us_K_fb",p.us_K)
    # save(filename,"N",p.basis.N,"Nz",p.basis.Nz,"a",p.basis.a,"L",p.basis.L,"v_dir",p.v_dir,"H_matrix",p.spec.H_matrix,"Es_K",p.Es_K,"us_K_fb",p.us_K)
    save_scfres(filename_root*".jld2", p.spec.scfres)
    px("Exported : V, u_K functions, for N=",p.basis.N,", Nz=",p.basis.Nz," at ",filename)
end

function import_u1_u2_V(N=20,Nz=120)
    p = Params(3)
    p.i_state = 4
    p.spec = InfoDFTK()
    p.basis.N = N; p.basis.Nz = Nz
    f_root = string(PATH_EXPORTS_MONOLAYER ,"N",p.basis.N,"_Nz",p.basis.Nz)
    f = f_root*".jld"

    p.basis.a = load(f,"a")
    p.basis.L = load(f,"L")
    # p.pseudo_inv = load(f,"pseudo_inv")
    # p.spec.H_matrix = load(f,"H_matrix")
    p.Es_K = load(f,"Es_K")
    p.us_K = load(f,"us_K_fb")

    p.spec.scfres = load_scfres(f_root*".jld2")

    init_cell_vectors(p.basis)
    init_cell_infinitesimals(p.basis)

    p.v_dir = load(f,"v_dir")
    # substract_by_far_value(p.v_dir,p) # don't do that !

    # Loads DFTK quantities
    p.spec.dftk_basis = basis_from_k(p.basis.K_red_3d,p)
    p.spec.K_kpt = p.spec.dftk_basis.kpoints[1]
    p.spec.Gplusk_vectors_cart = collect(DFTK.Gplusk_vectors_cart(p.spec.dftk_basis, p.spec.K_kpt))
    p.spec.Gvectors = collect(DFTK.G_vectors(p.spec.dftk_basis, p.spec.K_kpt))
    p.spec.Gvectors_inv = inverse_dict_from_array(Tuple.(p.spec.Gvectors)) # builds the inverse
    px("On imports : size of Fourier ball=",N_basis(p)," N=",N,", Nz=",Nz)
    p
end
