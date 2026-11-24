using LinearAlgebra, Statistics, StaticArrays, SparseArrays, Sunny
using Test

include(joinpath(@__DIR__, "../", "src/Phunny.jl"))
using .Phunny

# Silicon (diamond structure)
a = 5.431
lat = lattice_vectors(a, a, a, 90, 90, 90) # Sunny func
pos = [[0,0,0], [1/4,1/4,1/4]]
types = ["Si","Si"]
cryst = Crystal(lat, pos; types) # Sunny func

# Lookup masses and coherent scattering lengths
mass = mass_lookup(Symbol.(types)) 	     # Dict(String=>amu)
bcoh = bcoh_lookup(Symbol.(types))           # Vector per site (fm)

# Springs (toy Si-like constants)
kL = (i,j,rij)->175.0      # eV/Å^2
kT = (i,j,rij)->127.0

mdl = build_model(cryst; cutoff=2*a, use_sunny_radius=false,kL=kL, kT=kT)
Φ = assemble_force_constants!(mdl)

# 4D DSF grid
h = collect(range(-5.5, 5.5; length=101))
k = collect(range(-5.5, 5.5; length=101))
l = [0.0]
E = ω_grid(0.0, 200.0, 401)
T = 150.0
η = 3.0

S4 = onephonon_dsf_4d(mdl, Φ, h, k, l, E;
                             q_basis=:rlu, q_cell=:primitive,
                             cryst=cryst, T=T, σ=η,
                             mass_unit=:amu, bcoh=bcoh, threads=false)

# Unit tests
@testset "onephonon" begin
	println()
	@test size(S4) == (101, 101, 1, 401) 	# 4d
	@test all(isfinite, S4) 		# 4d
end

@testset "grid" begin
	println()
	@test E == 0.0:0.5:200.0
end

@testset "build_model" begin
	println()
	# the cutoff (2a) reaches neighbouring cells; the bond inside the home cell is still there
	b = only(filter(x -> x.i == 1 && x.j == 2 && iszero(x.R), mdl.bonds))
	@test b.r0 ≈ [1.35775, 1.35775, 1.35775]
	@test all(x -> norm(x.r0) <= 2a + 1e-9, mdl.bonds)
	@test any(x -> !iszero(x.R), mdl.bonds)
    @test mdl.lattice ≈ [5.431  0.0    0.0;
			 0.0    5.431  0.0;
 			 0.0    0.0    5.431]
    @test mdl.species == [:Si, :Si]
    @test mdl.N == 2
end;

@testset "force_constants" begin
	println()
	@test all(v -> size(v) == (3,3), values(Φ))
	@test haskey(Φ, (1, 1, SVector(0, 0, 0))) && haskey(Φ, (1, 2, SVector(0, 0, 0)))
	@test asr_residual(Φ, mdl.N) < 1e-8                  # translational invariance
end

@testset "mass_bcoh" begin
	println()
	@test mass ≈ [28.085, 28.085]
	@test bcoh ≈ [4.1491, 4.1491]
end

phonon_energies, phonon_vecs = phonons(mdl, Φ, @SVector[0.5,0.5,0.5]; q_basis=:cart, q_cell=:primitive)

include(joinpath(@__DIR__, "api_fixes.jl"))
include(joinpath(@__DIR__, "magnetoelastic.jl"))

@testset "Phunny Validation Suite" begin
    include(joinpath(@__DIR__, "Validation", "bcoh-test.jl"))
    include(joinpath(@__DIR__, "Validation", "U-tensor-test.jl"))
    include(joinpath(@__DIR__, "Validation", "DW-factor-test.jl"))
    include(joinpath(@__DIR__, "Validation", "onephonon-dsf-test.jl"))
end
print("Finished Tests!")
