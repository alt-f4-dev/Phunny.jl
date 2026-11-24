#-----------------------------------------#
# Zone-centre irreducible representations #
#-----------------------------------------#

"""
    displacement_representation(model, R, t; tol=1e-6) -> Matrix{Float64}

Matrix (3N × 3N) of the space-group operation {R | t} (`R`, `t` in fractional coordinates) on the Cartesian
displacements at Γ: the atom permutation it induces times the Cartesian rotation L R L⁻¹. Throws if the operation
does not map the structure onto itself or exchanges different species.
"""
function displacement_representation(model::Model, R::AbstractMatrix, t::AbstractVector; tol=1e-6)
    L = model.lattice; Rc = L * R / L; n = model.N
    G = zeros(3n, 3n)
    for i in 1:n
        r = R * model.fracpos[i] + t
        j = findfirst(f -> norm(mod.(r - f .+ 0.5, 1.0) .- 0.5) < tol, model.fracpos)
        (j === nothing || model.species[j] != model.species[i]) &&
            throw(ArgumentError("the operation does not map the structure onto itself"))
        G[3j-2:3j, 3i-2:3i] .= Rc
    end
    return G
end

"""
    gamma_irreps(model, D0, ops, table; skip=3, tol=1e-6) -> Vector{(E, irrep, χ)}

Zone-centre modes of the dynamical matrix `D0` = D(0) labelled by irreducible representation. `ops` is a list of
operations `(R, t)` (fractional), `table` maps each irrep name to its characters under `ops` (in the same order).
Degenerate modes are grouped and reported once; the characters χ of each group are matched to `table`
(`:unknown` if none matches within 1e-3). The lowest `skip` modes (the acoustic ones) are left out.
Eigenvectors of D are mass-weighted, which leaves characters unchanged since operations only exchange like atoms.
"""
function gamma_irreps(model::Model, D0::AbstractMatrix, ops, table::AbstractDict; skip::Integer=3, tol=1e-6)
    E, V = phonons(D0)
    Gs = [displacement_representation(model, R, t) for (R, t) in ops]
    out = Tuple{Float64,Symbol,Vector{Float64}}[]
    i = skip + 1
    while i <= length(E)
        j = i
        while j < length(E) && abs(E[j+1] - E[i]) < tol * max(1, abs(E[i])); j += 1; end
        U = V[:, i:j]
        χ = [real(tr(U' * G * U)) for G in Gs]
        name = something(findfirst(c -> length(c) == length(χ) && all(abs.(c .- χ) .< 1e-3), table), :unknown)
        push!(out, (E[i], name, χ))
        i = j + 1
    end
    return out
end
