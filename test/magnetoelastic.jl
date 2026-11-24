# Fitting, zone-centre symmetry and magnon–phonon hybridization on a two-atom cubic toy (CsCl structure).
using LinearAlgebra, Random, StaticArrays, Statistics, Test

const TOY = (lattice=Matrix(3.0I, 3, 3), positions=[[0.0, 0, 0], [0.5, 0.5, 0.5]], species=[:A, :B])
toy_model() = build_model(TOY; mass=Dict(:A => 10.0, :B => 20.0), cutoff=3.0, use_sunny_radius=false)
const TOY_GROUPS = Dict((:A, :B) => :AB, (:A, :A) => :AA, (:B, :B) => :BB)
const TOY_K = [0.3, 0.05, 1.2, 0.2, 0.4, 0.1]

@testset "Fitting force constants" begin
    fb = ForceConstantBasis(toy_model(), Dict(1 => :A, 2 => :B), TOY_GROUPS)
    Dp = dynamical_matrices(fb, SVector(0.3, -0.2, 0.5))
    @test dynamical_matrix(Dp, TOY_K) ≈ dynamical_matrix(fb, TOY_K, SVector(0.3, -0.2, 0.5))
    θ1 = debye_temperature(fb, TOY_K)
    @test isfinite(θ1)
    @test debye_temperature(fb, 4 .* TOY_K) ≈ 2θ1 rtol = 1e-8     # θ ∝ √k
    qs = [SVector(0.4, 0.1, 0.0), SVector(0.2, 0.3, 0.5), SVector(0.5, 0.5, 0.1)]
    obs(k) = (E = reduce(vcat, [phonons(dynamical_matrix(fb, k, q))[1] for q in qs]); minimum(E) > 0 ? E : nothing)
    data = obs(TOY_K); σ = 0.01 .* data
    fit = fit_force_constants(obs, data, σ; lower=fill(1e-3, 6), upper=fill(10.0, 6), starts=60,
                              rng=MersenneTwister(1))
    @test fit.k ≈ TOY_K rtol = 1e-4
    @test fit.χ2ν < 1e-6 && fit.model ≈ data
    s, acc = sample_force_constants(obs, data, σ, fit.θ, 50; lower=fill(1e-3, 6), upper=fill(10.0, 6), step=0.005,
                                    rng=MersenneTwister(2))
    @test size(s) == (6, 50) && 0.05 < acc < 0.95
    @test exp.(vec(mean(s; dims=2))) ≈ TOY_K rtol = 0.05
end

@testset "Zone-centre irreps" begin
    mdl = toy_model(); fb = ForceConstantBasis(mdl, Dict(1 => :A, 2 => :B), TOY_GROUPS)
    inv = (Matrix(-1.0I, 3, 3), zeros(3)); c4 = ([0.0 -1 0; 1 0 0; 0 0 1], zeros(3))
    @test displacement_representation(mdl, inv...)' * displacement_representation(mdl, inv...) ≈ I
    @test_throws ArgumentError displacement_representation(mdl, Matrix(1.0I, 3, 3), [0.5, 0.5, 0.5])  # A ↔ B
    g = gamma_irreps(mdl, dynamical_matrix(fb, TOY_K, SVector(0.0, 0, 0)), [inv, c4], Dict(:T1u => [-3.0, 1.0]))
    @test length(g) == 1 && g[1][2] === :T1u                # the optical triplet is polar
    @test gamma_irreps(mdl, dynamical_matrix(fb, TOY_K, SVector(0.0, 0, 0)), [inv], Dict(:T1g => [3.0]))[1][2] === :unknown
end

@testset "Magnon–phonon hybridization" begin
    # colpa: one boson, ω = √(A² − B²), para-unitary
    E1, T1 = colpa([2.0 1.0; 1.0 2.0])
    @test E1[1] ≈ sqrt(3)
    @test T1' * Diagonal([1, -1]) * T1 ≈ I
    mdl = toy_model(); fb = ForceConstantBasis(mdl, Dict(1 => :A, 2 => :B), TOY_GROUPS)
    env = LocalEnvironment(mdl, 1; rmax=2.7, include=κ -> mdl.species[κ] === :B)
    @test length(env.neighbors) == 8 && all(n -> norm(n[2]) ≈ 1.5sqrt(3), env.neighbors)
    G0 = randn(MersenneTwister(3), 3, 3)
    @test G0 * reduce(hcat, last.(env.neighbors)) * env.P ≈ G0           # affine fields are recovered exactly
    # a long acoustic wave has the gradient i x₀ ê qᵀ: rank one, annihilating vectors ⊥ q
    q = SVector(0.004, 0.0, 0.0); E, e = phonons(dynamical_matrix(fb, TOY_K, q))
    G = displacement_gradient(mdl, env, q, e[:, 1], E[1])
    @test opnorm(G) ≈ mode_amplitude(30.0, E[1]) * norm(q) rtol = 1e-3
    @test norm(G * [0, 1, 0]) + norm(G * [0, 0, 1]) < 1e-6 * opnorm(G)
    # the single-ion coupling: linear in B, and a pure rotation enters through D alone
    f = single_ion_magnetoelastic(B=2.0, D=0.5)
    ω = [0 0 0.1; 0 0 0; -0.1 0 0]; ε = [0 0 0.2; 0 0 0; 0.2 0 0]
    @test f(zeros(3, 3), ω, 2.0) == (-3 * 0.5 * 0.1, 0.0)
    @test f(ε, zeros(3, 3), 2.0)[1] ≈ 2 * single_ion_magnetoelastic(B=1.0)(ε, zeros(3, 3), 2.0)[1]
    @test single_ion_magnetoelastic(B=0.0, D=0.5, rotational=false)(ε, ω, 2.0) == (0.0, 0.0)
    # one magnetic site tuned onto a phonon at a generic q
    site = MagneticSite(env, SMatrix{3,3}(1.0I), 2.0)
    q = SVector(0.31, 0.17, 0.08); ph = phonons(dynamical_matrix(fb, TOY_K, q))
    ν = 4; a = ph[1][ν]; Mmag = [a 0.0; 0.0 a]
    M0 = magnon_phonon_bdg(Mmag, mdl, [site], q, ph)
    @test hybrid_modes(M0, 1)[1] ≈ sort([ph[1]; a])
    weak = single_ion_magnetoelastic(B=0.02)
    M = magnon_phonon_bdg(Mmag, mdl, [site], q, ph; coupling=weak)
    Eh, mag = hybrid_modes(M, 1)
    @test sum(mag) ≈ 1
    V = first_order_coupling(M, M0, 1)
    @test size(V) == (1, 6) && abs(V[1, ν]) > 0
    near = sort(Eh[sortperm(abs.(Eh .- a))[1:2]])
    @test near[2] - near[1] ≈ 2abs(V[1, ν]) rtol = 1e-2
    @test first_order_coupling(magnon_phonon_bdg(Mmag, mdl, [site], q, ph; coupling=single_ion_magnetoelastic(B=0.04)), M0, 1) ≈ 2V
    _, T = colpa(magnon_phonon_bdg(Mmag, mdl, [site], q, ph; coupling=single_ion_magnetoelastic(B=5.0)))
    @test T' * Diagonal([ones(7); -ones(7)]) * T ≈ I
    @test_throws DomainError magnon_phonon_bdg(Mmag, mdl, [site], SVector(0.0, 0, 0),
                                               phonons(dynamical_matrix(fb, TOY_K, SVector(0.0, 0, 0))))
end
