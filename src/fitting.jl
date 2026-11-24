#--------------------------------------------------------------#
# Fitting force constants to measured data, and their posterior #
#--------------------------------------------------------------#

"""
    dynamical_matrix(Dp::AbstractArray{<:Complex,3}, k) -> Hermitian

Σₚ kₚ Dₚ for basis matrices cached with [`dynamical_matrices`](@ref).
"""
function dynamical_matrix(Dp::AbstractArray{<:Complex,3}, k::AbstractVector)
    size(Dp, 3) == length(k) || throw(DimensionMismatch("need one force constant per basis matrix"))
    D = zeros(ComplexF64, size(Dp, 1), size(Dp, 2))
    @inbounds for p in eachindex(k)
        D .+= k[p] .* @view Dp[:, :, p]
    end
    return Hermitian(D)
end

# Fibonacci directions on the upper hemisphere (sound velocities are even in q).
_hemisphere(n) = [(φ = π * (3 - √5) * i; z = 1 - (i + 0.5) / n; r = sqrt(1 - z^2); SVector(r * cos(φ), r * sin(φ), z))
                  for i in 0:n-1]

"""
    debye_temperature(fb, k; directions=40, dq=2e-3) -> K

Elastic Debye temperature θ = (ħ/k_B) v_m (6π² n)^{1/3} of the force constants `k`, with the mean sound velocity
v_m⁻³ = ⟨Σₛ vₛ⁻³⟩/3 over `directions` directions (finite difference at |q| = `dq` Å⁻¹) and n the atomic density.
`NaN` if an acoustic branch is unstable.
"""
function debye_temperature(fb::ForceConstantBasis, k::AbstractVector; directions::Integer=40, dq::Real=2e-3)
    inv3 = 0.0; dirs = _hemisphere(directions)
    for d in dirs
        E, _ = phonons(dynamical_matrix(fb, k, SVector{3,Float64}(dq .* d)))
        any(<=(0), E[1:3]) && return NaN
        inv3 += sum((151.927 .* E[1:3] ./ dq) .^ -3) / 3          # m/s from meV and Å⁻¹
    end
    vm = (inv3 / length(dirs))^(-1 / 3)
    n = fb.model.N / (abs(det(fb.model.lattice)) * 1e-30)
    return 7.63823e-12 * vm * cbrt(6π^2 * n)
end

"""
    fit_force_constants(observables, data, σ; lower, upper, starts=300, best=10, iters=80, rng)
        -> (; θ, k, cov, χ2ν, residuals, model, data)

Weighted least-squares fit of positive parameters k (typically basis force constants) to `data ± σ`.
`observables(k)` returns the model values in the order of `data`, or `nothing` if `k` is unphysical (e.g. an
unstable lattice). The fit runs in θ = log k inside the box [log `lower`, log `upper`]: `starts` uniform random
points are screened, and Levenberg–Marquardt refines the `best` of them. `cov` is the covariance of θ, scaled by
χ²ν when it exceeds 1 (pseudo-inverse, so unconstrained directions do not fail).
"""
function fit_force_constants(observables, data::AbstractVector, σ::AbstractVector; lower::AbstractVector,
                             upper::AbstractVector, starts::Integer=300, best::Integer=10, iters::Integer=80,
                             rng=Random.default_rng())
    lo, hi = log.(lower), log.(upper); np = length(lo)
    function resid(θ)
        o = observables(exp.(θ))
        return o === nothing ? fill(1e3, length(data)) : (o .- data) ./ σ
    end
    function lm(θ)
        r = resid(θ); χ2 = sum(abs2, r); λ = 1e-2; J = zeros(length(r), np)
        for _ in 1:iters
            for j in 1:np
                δ = zeros(np); δ[j] = 1e-4
                J[:, j] .= (resid(θ .+ δ) .- r) ./ 1e-4
            end
            A = J' * J; g = J' * r; improved = false
            for _ in 1:10
                trial = clamp.(θ .- (A + λ * Diagonal(diag(A) .+ 1e-12)) \ g, lo, hi)
                rt = resid(trial); χt = sum(abs2, rt)
                if χt < χ2
                    improved = χ2 - χt >= 1e-8 * χ2
                    θ, r, χ2 = trial, rt, χt; λ = max(λ / 5, 1e-10)
                    break
                end
                λ *= 8
            end
            improved || break
        end
        return θ, r, χ2, J
    end
    cands = [lo .+ (hi .- lo) .* rand(rng, np) for _ in 1:starts]
    order = sortperm([sum(abs2, resid(c)) for c in cands])
    fits = [lm(c) for c in cands[order[1:min(best, end)]]]
    θ, r, χ2, J = fits[argmin([f[3] for f in fits])]
    ν = max(length(data) - np, 1)
    return (; θ, k=exp.(θ), cov=pinv(J' * J) .* max(χ2 / ν, 1.0), χ2ν=χ2 / ν, residuals=r,
            model=observables(exp.(θ)), data=collect(data))
end

"""
    sample_force_constants(observables, data, σ, θ0, n; lower, upper, step=0.15, thin=10, burn=300, rng)
        -> (samples, acceptance)

Random-walk Metropolis over θ = log k with likelihood exp(−χ²/2), uniform prior in the box of
[`fit_force_constants`](@ref), started at `θ0` (e.g. the fit). Returns `n` thinned samples (columns of θ) of the
posterior and the acceptance rate. Unphysical proposals (`observables` returns `nothing`) are rejected.
"""
function sample_force_constants(observables, data::AbstractVector, σ::AbstractVector, θ0::AbstractVector, n::Integer;
                                lower::AbstractVector, upper::AbstractVector, step::Real=0.15, thin::Integer=10,
                                burn::Integer=300, rng=Random.default_rng())
    lo, hi = log.(lower), log.(upper)
    χ2(θ) = (o = observables(exp.(θ)); o === nothing ? Inf : sum(abs2, (o .- data) ./ σ))
    θ = collect(float.(θ0)); c = χ2(θ)
    out = Matrix{Float64}(undef, length(θ), n); kept = accepted = total = 0
    while kept < n
        trial = θ .+ step .* randn(rng, length(θ))
        if all(lo .<= trial .<= hi)
            ct = χ2(trial)
            log(rand(rng)) < (c - ct) / 2 && ((θ, c) = (trial, ct); accepted += 1)
        end
        total += 1
        total > burn && (total - burn) % thin == 0 && (kept += 1; out[:, kept] .= θ)
    end
    return out, accepted / total
end
