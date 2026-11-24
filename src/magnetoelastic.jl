#---------------------------------------------------------------#
# Magnon–phonon hybridization of collinear magnets               #
#---------------------------------------------------------------#
#
# Conventions. Phonon eigenvectors e(q) are those of `dynamical_matrix` (atom-position phases e^{iq·r}); the
# displacement of atom κ is u_κ = x₀(M_κ, E) e_κ(q) per unit of (c_q + c†₋q), x₀ = `mode_amplitude`. Magnon
# blocks must use the same atom-position phases. The BdG basis is Ψ = (a_q, c_q, a†₋q, c†₋q), with the magnon
# bosons a first, and H = ½ Ψ† M Ψ, M = [A(q) B(q); B(q)† A(−q)*]. Phonons at −q use e(−q) = e(q)*.
#
# Each magnetic site carries a local frame F (columns x̂′, ŷ′, ẑ′; ẑ′ along the ordered moment). To linear order
# in the Holstein–Primakoff bosons S^{x′} = √(S/2)(a + a†), S^{y′} = −i√(S/2)(a − a†), so a coupling written as
# h_x S^{x′} + h_y S^{y′} contributes √(S/2)(h_x + i h_y) to the coefficient of a†.

"""
    colpa(M) -> (E, T)

Bosonic Bogoliubov diagonalization (Colpa, Physica A 93, 327 (1978)) of a positive-definite 2n × 2n BdG matrix.
Returns the n positive energies (ascending) and the corresponding 2n × n columns of the para-unitary
transformation, normalized so that `T' * σ₃ * T = I` with σ₃ = diag(1ₙ, −1ₙ). A `PosDefException` means the
quadratic Hamiltonian is not bounded below: the assumed ordered state is unstable.
"""
function colpa(M::AbstractMatrix)
    n = size(M, 1) ÷ 2
    K = cholesky(Hermitian(Matrix(M))).U
    vals, U = eigen!(Hermitian(K * Diagonal([ones(n); -ones(n)]) * K'))
    pos = n+1:2n
    return vals[pos], K \ (U[:, pos] .* sqrt.(vals[pos])')
end

"""
    LocalEnvironment(model, center; rmax, include=κ -> true)

The neighbours of atom `center` within `rmax` (Å, all periodic images) that pass `include(κ)`, with the projector
P = Δᵀ(ΔΔᵀ)⁻¹ of their bond vectors Δ (3 × n). For relative displacements U (3 × n) the least-squares
displacement gradient is G = U P, whose symmetric and antisymmetric parts are the local strain and rotation.
At least three non-coplanar neighbours are required.
"""
struct LocalEnvironment
    center::Int
    neighbors::Vector{Tuple{Int,SVector{3,Float64}}}
    P::Matrix{Float64}
end

function LocalEnvironment(model::Model, center::Integer; rmax::Real, include=κ -> true)
    L = model.lattice; f0 = model.fracpos[center]
    reach = ceil.(Int, rmax .* [norm(cross(L[:, mod1(i + 1, 3)], L[:, mod1(i + 2, 3)])) / abs(det(L)) for i in 1:3]) .+ 1
    neigh = Tuple{Int,SVector{3,Float64}}[]
    for κ in 1:model.N, R in Iterators.product((-r:r for r in reach)...)
        δ = L * (model.fracpos[κ] + SVector{3,Float64}(R) - f0)
        0 < norm(δ) <= rmax && include(κ) && push!(neigh, (κ, δ))
    end
    Δ = reduce(hcat, last.(neigh); init=zeros(3, 0))
    rank(Δ) == 3 || throw(ArgumentError("atom $center needs three non-coplanar neighbours within $rmax Å"))
    return LocalEnvironment(center, neigh, Δ' / (Δ * Δ'))
end

"""
    displacement_gradient(model, env, q, e, E) -> Matrix{ComplexF64}

Least-squares displacement gradient G = ∂u/∂x (3 × 3) of the environment `env` in the phonon mode of energy `E`
(meV) and eigenvector `e` (length 3N) at the Cartesian wave vector `q`, per unit of (c + c†). Phases are taken
relative to the central atom.
"""
function displacement_gradient(model::Model, env::LocalEnvironment, q, e::AbstractVector, E::Real)
    c = env.center
    u0 = mode_amplitude(model.mass[c], E) .* SVector(e[3c-2], e[3c-1], e[3c])
    U = zeros(ComplexF64, 3, length(env.neighbors))
    for (n, (κ, δ)) in enumerate(env.neighbors)
        U[:, n] .= mode_amplitude(model.mass[κ], E) .* cis(dot(q, δ)) .* SVector(e[3κ-2], e[3κ-1], e[3κ]) .- u0
    end
    return U * env.P
end

"""
    MagneticSite(env, frame, S)

A magnetic atom (the centre of `env`) with spin `S` and local frame `frame` (3 × 3, columns x̂′, ŷ′, ẑ′ in
Cartesian coordinates, ẑ′ along the ordered moment). Sites are listed in the order of the magnon bosons.
"""
struct MagneticSite
    env::LocalEnvironment
    frame::SMatrix{3,3,Float64,9}
    S::Float64
end

"""
    single_ion_magnetoelastic(; B, D=0.0, rotational=true)

Coupling of a uniaxial single ion (−D S_z′² about the moment) to the strain ε and rotation ω of its environment,
H = Σ_{α=x′,y′} (B ε_{αz′} − D ω_{αz′}) O_{αz′}, with O_{αz′} = S_α S_z′ + S_z′ S_α. The rotational part is the
one fixed by rotational invariance (Melcher, PRL 28, 165 (1972)); `rotational=false` drops it. To linear order
O_{αz′} → (2S − 1) S_α. Returns `(ε′, ω′, S) -> (h_x, h_y)`, the fields conjugate to S^{x′}, S^{y′}, with ε′, ω′ in
the local frame. Any function of this form can be passed to [`magnon_phonon_bdg`](@ref) instead.
"""
function single_ion_magnetoelastic(; B::Real, D::Real=0.0, rotational::Bool=true)
    Dr = rotational ? D : 0.0
    return (ε, ω, S) -> ((2S - 1) * (B * ε[1, 3] - Dr * ω[1, 3]), (2S - 1) * (B * ε[2, 3] - Dr * ω[2, 3]))
end

# √(S/2)(h_x + i h_y) of one site in one phonon mode at q (eigenvector e).
function _site_coupling(model, site, coupling, q, e, E)
    G = displacement_gradient(model, site.env, q, e, E)
    F = site.frame
    ε = F' * ((G + transpose(G)) / 2) * F
    ω = F' * ((G - transpose(G)) / 2) * F
    hx, hy = coupling(ε, ω, site.S)
    return sqrt(site.S / 2) * (hx + im * hy)
end

"""
    magnon_phonon_bdg(Mmag, model, sites, q, (E, e); coupling=nothing, Emin=1e-3) -> Matrix{ComplexF64}

BdG matrix of the coupled magnons and phonons at the Cartesian wave vector `q`. `Mmag` is the 2m × 2m magnon BdG
matrix [A(q) B(q); B(q)† A(−q)*] (m = `length(sites)`, atom-position phases), `(E, e)` the phonon energies (all
> `Emin` meV: the zero-point amplitude diverges for the acoustic modes at Γ) and eigenvectors from [`phonons`](@ref). `coupling(ε′, ω′, S) -> (h_x, h_y)` is applied at every
site (see [`single_ion_magnetoelastic`](@ref)); with `coupling = nothing` the blocks are uncoupled.
"""
function magnon_phonon_bdg(Mmag::AbstractMatrix, model::Model, sites::AbstractVector{MagneticSite}, q, ph;
                           coupling=nothing, Emin=1e-3)
    E, e = ph
    minimum(E) > Emin || throw(DomainError(minimum(E), "phonon energies must exceed Emin = $Emin meV (avoid Γ)"))
    m = length(sites); size(Mmag) == (2m, 2m) || throw(DimensionMismatch("magnon BdG must be $(2m) × $(2m)"))
    nph = length(E); n = m + nph; q = SVector{3,Float64}(q)
    A = zeros(ComplexF64, n, n); Am = zeros(ComplexF64, n, n); Bm = zeros(ComplexF64, n, n)
    A[1:m, 1:m] .= Mmag[1:m, 1:m]; Bm[1:m, 1:m] .= Mmag[1:m, m+1:2m]; Am[1:m, 1:m] .= conj.(Mmag[m+1:2m, m+1:2m])
    for ν in 1:nph
        A[m+ν, m+ν] = Am[m+ν, m+ν] = E[ν]
    end
    if coupling !== nothing
        for (i, s) in enumerate(sites), ν in 1:nph
            ev = @view e[:, ν]
            G = _site_coupling(model, s, coupling, q, ev, E[ν])
            Gm = _site_coupling(model, s, coupling, -q, conj.(ev), E[ν])
            A[i, m+ν] = G; A[m+ν, i] = conj(G)
            Am[i, m+ν] = Gm; Am[m+ν, i] = conj(Gm)
            Bm[i, m+ν] = G; Bm[m+ν, i] = Gm
        end
    end
    return [A Bm; Bm' conj(Am)]
end

"""
    hybrid_modes(M, m) -> (E, magnon_fraction)

Energies (ascending) of the BdG matrix `M` and, for each mode, its para-norm weight on the first `m` bosons
(the magnons): Σᵢ≤ₘ |Tᵢⱼ|² − |T_{n+i,j}|². The fractions sum to m.
"""
function hybrid_modes(M::AbstractMatrix, m::Integer)
    E, T = colpa(M); n = size(M, 1) ÷ 2
    return E, [sum(abs2, @view T[1:m, j]) - sum(abs2, @view T[n+1:n+m, j]) for j in 1:n]
end

"""
    first_order_coupling(M, M0, m) -> Matrix{ComplexF64}

Matrix elements t_j† (M − M₀) t_ν between the normal modes of the uncoupled magnons (rows, ascending energy) and
the phonons (columns, in the order of `M0`), where M₀ is the uncoupled matrix (`coupling = nothing`). Near a
crossing of a magnon and a phonon detuned by δ, the exact splitting is √(δ² + 4|V|²) to leading order.
"""
function first_order_coupling(M::AbstractMatrix, M0::AbstractMatrix, m::Integer)
    n = size(M0, 1) ÷ 2
    _, Tm = colpa(M0[[1:m; n+1:n+m], [1:m; n+1:n+m]])
    t = zeros(ComplexF64, 2n, m); t[1:m, :] .= Tm[1:m, :]; t[n+1:n+m, :] .= Tm[m+1:2m, :]
    return (t' * (M - M0))[:, m+1:n]
end
