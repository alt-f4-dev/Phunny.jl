# Neighbour search across cells, bond deduplication, signed energies, the linear force-constant basis, amplitudes.
using LinearAlgebra, StaticArrays, Test

@testset "Neighbours across cell boundaries" begin
    L = SMatrix{3,3,Float64,9}(3.0I)
    b = neighbor_bonds_cutoff(L, [SVector(0.0, 0.0, 0.0)]; cutoff=3.1)
    @test length(b) == 3                                   # ±x, ±y, ±z, one bond per pair
    @test sort([Tuple(x.R) for x in b]) == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    # face diagonals with mixed-sign offsets survive deduplication: 3 axes + 6 face diagonals
    b2 = neighbor_bonds_cutoff(L, [SVector(0.0, 0.0, 0.0)]; cutoff=4.3)
    @test length(b2) == 9 && any(x -> Tuple(x.R) == (1, -1, 0), b2)
    # a long cutoff in a small cell reaches beyond the nearest image
    @test length(neighbor_bonds_cutoff(L, [SVector(0.0, 0.0, 0.0)]; cutoff=6.1)) == 3 + 6 + 4 + 3
end

@testset "Signed energies and the force-constant basis" begin
    cell = (lattice=Matrix(3.0I, 3, 3), positions=[[0.0, 0, 0], [0.5, 0.5, 0.5]], species=[:A, :B])
    mdl = build_model(cell; mass=Dict(:A => 10.0, :B => 20.0), cutoff=3.0, use_sunny_radius=false)
    atom = Dict(1 => :A, 2 => :B)
    fb = ForceConstantBasis(mdl, atom, Dict((:A, :B) => :AB, (:A, :A) => :AA, (:B, :B) => :BB))
    @test fb.params == [(:AA, :L), (:AA, :T), (:AB, :L), (:AB, :T), (:BB, :L), (:BB, :T)]
    k = [0.3, 0.05, 1.2, 0.2, 0.4, 0.1]
    assign_force_constants!(mdl, atom, Dict((:B, :A) => (1.2, 0.2), (:A, :A) => (0.3, 0.05), (:B, :B) => (0.4, 0.1)))
    q = SVector(0.3, -0.2, 0.5)
    @test dynamical_matrix(fb, k, q) ≈ dynamical_matrix(mdl, assemble_force_constants!(mdl), q)   # either pair order
    E, _ = phonons(Matrix(dynamical_matrix(fb, k, q)))
    @test all(>(0), E)
    Eu, _ = phonons(Matrix(dynamical_matrix(fb, [0.3, 0.05, -1.2, 0.2, 0.4, 0.1], q)))
    @test minimum(Eu) < 0                                  # an unstable mode is reported, not hidden
    @test minimum(phonons(Matrix(dynamical_matrix(fb, [0.3, 0.05, -1.2, 0.2, 0.4, 0.1], q)); signed=false)[1]) == 0
    @test mode_amplitude(1.0, HBAR2_AMU_A2 / 2) ≈ 1.0
    # groups by bond length: same matrices as by label pairs when the lengths separate the same bonds
    byd = ForceConstantBasis(mdl, b -> atom[b.i] == atom[b.j] ? (atom[b.i] === :A ? :AA : :BB) : :AB)
    @test byd.params == fb.params && dynamical_matrix(byd, k, q) ≈ dynamical_matrix(fb, k, q)
    @test_throws DimensionMismatch dynamical_matrix(fb, k[1:3], q)
    # bond-bending terms are linear too
    mb = build_model(cell; mass=Dict(:A => 10.0, :B => 20.0), cutoff=3.0, use_sunny_radius=false, β=0.7)
    @test_throws ArgumentError ForceConstantBasis(mb, b -> :all)
    fa = ForceConstantBasis(mb, b -> :all; angle_group=a -> :bend)
    @test fa.params[end] == (:bend, :β)
    assign_force_constants!(mb, atom, Dict((:A, :B) => (1.0, 0.1), (:A, :A) => (1.0, 0.1), (:B, :B) => (1.0, 0.1)))
    @test dynamical_matrix(fa, [1.0, 0.1, 0.7], q) ≈ dynamical_matrix(mb, assemble_force_constants!(mb), q)
end
