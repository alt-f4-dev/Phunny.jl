
#These are to be added into Physical.jl (top-level) once testing is finished
#
#export DisorderSpec, DisorderModel, DisorderEnsemble
#export make_supercell, dope, vacate, distort
#export build_disordered_model, sample_disordered_model

#-----------------------------------------------#
#               Core abstractions               #
#-----------------------------------------------#

"""
    DisorderSpec

Specification of an allowed disorder family.

This is a declarative object. It does not represent one realized crystal.
"""
Base.@kwdef struct DisorderSpec
    dims::NTuple{3,Int} = (1,1,1)
    
    #--------------------------------------------------------------#
    # substitution rules                                           #
    #                                                              #
    # Example:                                                     #
    #   Dict(:Cu => (to=:Zn, fraction=0.10, allowed_sites=nothing))#
    #--------------------------------------------------------------#
    substitutions::Dict{Symbol,Any} = Dict{Symbol,Any}()

    #--------------------------------------------------------------#
    # vacancy rules                                                #
    #                                                              #
    # Example:                                                     #
    #   Dict(:O => (fraction=0.05, allowed_sites=nothing))         #
    #--------------------------------------------------------------#
    vacancies::Dict{Symbol,Any} = Dict{Symbol,Any}()

    #--------------------------------------------------------------#
    # distortion rule on fractional positions or cartesian         #
    # positions should accept (lattice, fracpos, species, rng)     #
    # and return new fracpos                                       #
    #--------------------------------------------------------------#
    displacement_rule::Union{Nothing,Function} = nothing

    #reproducibility
    seed::Union{Nothing,Int} = nothing

    #optional metadata for force-constant provenance
    label::String = "disorder_spec"
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end


"""
    DisorderModel

One explicit disorder realization and its built Phunny model.
This is your realized configuration object.
"""
Base.@kwdef mutable struct DisorderModel
    parent_crystal
    realized_crystal
    spec::DisorderSpec

    #structure data
    dims::NTuple{3,Int} = (1,1,1)
    fracpos::Vector = Any[]
    species::Vector{Symbol} = Symbol[]
    site_map::Vector{Int} = Int[]          # map back to parent-site labels if useful
    substitutions::Vector{NamedTuple} = NamedTuple[]
    vacancies::Vector{NamedTuple} = NamedTuple[]
    seed::Union{Nothing,Int} = nothing

    #built forward model + results cache
    model = nothing
    Φ = nothing

    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end


"""
    DisorderEnsemble

Collection of many disorder realizations plus ensemble statistics.
"""
Base.@kwdef mutable struct DisorderEnsemble
    parent_crystal
    spec::DisorderSpec
    nsamples::Int = 0
    samples::Vector{DisorderModel} = DisorderModel[]

    # optional ensemble-level summaries
    observables::Vector{Any} = Any[]
    mean_observable = nothing
    var_observable = nothing

    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
end


#-------------------------------------------------------#
#                       Utilities                       #
#-------------------------------------------------------#

# Extract crystal fields in the same broad spirit as Phunny._to_phunny_spec
function _extract_crystal(cryst)
    if hasproperty(cryst, :latvecs) && hasproperty(cryst, :positions) && hasproperty(cryst, :types)
        L  = SMatrix{3,3,Float64,9}(Matrix(getproperty(cryst, :latvecs)))
        fp = [SVector{3,Float64}(Tuple(p)) for p in getproperty(cryst, :positions)]
        sp = Symbol.(getproperty(cryst, :types))
        return (lattice=L, fracpos=fp, species=sp)
    elseif hasproperty(cryst, :lattice) && hasproperty(cryst, :positions) && hasproperty(cryst, :species)
        L  = SMatrix{3,3,Float64,9}(Matrix(getproperty(cryst, :lattice)))
        fp = [SVector{3,Float64}(Tuple(p)) for p in getproperty(cryst, :positions)]
        sp = Symbol.(getproperty(cryst, :species))
        return (lattice=L, fracpos=fp, species=sp)
    else
        error("Unsupported crystal object interface.")
    end
end


# Replace this with your preferred crystal constructor.
# For Sunny users this can be `Crystal(L, fracpos; types=String.(species))`.
function _rebuild_crystal(lattice, fracpos, species)
    return (; latvecs=Matrix(lattice), positions=fracpos, types=String.(species))
end


#---------------------------------------------------------------#
#               1. Supercell Construction                       #
#---------------------------------------------------------------#

"""
    make_supercell(cryst, dims)

Build an explicit repeated supercell crystal from `cryst`.
Returns a new crystal-like object.
"""
function make_supercell(cryst, dims::NTuple{3,Int})
    nx, ny, nz = dims
    spec = _extract_crystal(cryst)
    L0, f0, s0 = spec.lattice, spec.fracpos, spec.species

    # new lattice columns scaled by dims
    S = SMatrix{3,3,Float64,9}(nx,0,0, 0,ny,0, 0,0,nz)
    Ls = L0 * S

    fs = SVector{3,Float64}[]
    ss = Symbol[]
    site_map = Int[]

    for iz in 0:nz-1, iy in 0:ny-1, ix in 0:nx-1
        T = SVector{3,Float64}(ix, iy, iz)
        for (n, f) in enumerate(f0)
            push!(fs, (f + T) ./ SVector{3,Float64}(nx, ny, nz))
            push!(ss, s0[n])
            push!(site_map, n)
        end
    end

    sc = _rebuild_crystal(Ls, fs, ss)
    return sc, site_map
end


#---------------------------------------------------------------#
#       2. Disorder Operations on Explicit Crystals             #
#---------------------------------------------------------------#

"""
    dope(cryst, dopant; rng=Random.default_rng())

Substitute selected sites according to `dopant`.

Expected `dopant` format examples:
    (from=:Cu, to=:Zn, fraction=0.1, allowed_sites=nothing)
    (from=:O,  to=:F,  fraction=0.2, allowed_sites=[1,4,7,...])

Returns:
    new_crystal, substitution_log
"""
function dope(cryst, dopant; rng=Random.default_rng())
    spec = _extract_crystal(cryst)
    L, fpos, species = spec.lattice, copy(spec.fracpos), copy(spec.species)

    from = dopant.from
    to = dopant.to
    frac = dopant.fraction
    allowed_sites = get(dopant, :allowed_sites, nothing)

    candidate_sites = findall(==(from), species)
    if allowed_sites !== nothing
        aset = Set(allowed_sites)
        candidate_sites = [i for i in candidate_sites if i in aset]
    end

    nsub = round(Int, frac * length(candidate_sites))
    chosen = nsub == 0 ? Int[] : Random.randperm(rng, length(candidate_sites))[1:nsub]
    chosen_sites = isempty(chosen) ? Int[] : candidate_sites[chosen]

    subs = NamedTuple[]
    for i in chosen_sites
        old = species[i]
        species[i] = to
        push!(subs, (site=i, from=old, to=to))
    end

    newcryst = _rebuild_crystal(L, fpos, species)
    return newcryst, subs
end


"""
    vacate(cryst, vacancy_rule; rng=Random.default_rng())

Remove selected sites.

Expected `vacancy_rule` format examples:
    (species=:O, fraction=0.05, allowed_sites=nothing)

Returns:
    new_crystal, vacancy_log
"""
function vacate(cryst, vacancy_rule; rng=Random.default_rng())
    spec = _extract_crystal(cryst)
    L, fpos, species = spec.lattice, copy(spec.fracpos), copy(spec.species)

    target = vacancy_rule.species
    frac = vacancy_rule.fraction
    allowed_sites = get(vacancy_rule, :allowed_sites, nothing)

    candidate_sites = findall(==(target), species)
    if allowed_sites !== nothing
        aset = Set(allowed_sites)
        candidate_sites = [i for i in candidate_sites if i in aset]
    end

    nvac = round(Int, frac * length(candidate_sites))
    chosen = nvac == 0 ? Int[] : Random.randperm(rng, length(candidate_sites))[1:nvac]
    chosen_sites = isempty(chosen) ? Int[] : sort(candidate_sites[chosen])

    vacs = NamedTuple[]
    keep = trues(length(species))
    for i in chosen_sites
        keep[i] = false
        push!(vacs, (site=i, species=species[i]))
    end

    newf = [fpos[i] for i in eachindex(fpos) if keep[i]]
    news = [species[i] for i in eachindex(species) if keep[i]]

    newcryst = _rebuild_crystal(L, newf, news)
    return newcryst, vacs
end


"""
    distort(cryst, displacement_rule; rng=Random.default_rng())

Apply an explicit distortion rule to atomic positions.

`displacement_rule` should return new fractional positions, or cartesian
displacements that you convert inside the rule.
"""
function distort(cryst, displacement_rule; rng=Random.default_rng())
    spec = _extract_crystal(cryst)
    L, fpos, species = spec.lattice, copy(spec.fracpos), copy(spec.species)

    newf = displacement_rule(L, fpos, species, rng)
    newcryst = _rebuild_crystal(L, newf, species)
    return newcryst
end


#-----------------------------------------------------------------------#
#                       3. Model Sample Builder                         #
#-----------------------------------------------------------------------#

"""
    build_disordered_model(cryst, spec; kwargs...)

Construct one explicit disorder realization, then feed it through
Phunny.build_model and optionally assemble Φ.

Keyword arguments are forwarded to Phunny.build_model.
"""
function build_disordered_model(cryst, spec::DisorderSpec;
                                rng = spec.seed === nothing ? Random.default_rng() : MersenneTwister(spec.seed),
                                assemble_Φ::Bool = true,
                                kL=0.0, kT=0.0, β=0.0,
                                mass=:lookup,
                                isotopes_by_site=nothing,
                                isotopes_by_species=nothing,
                                neighbors_sunny=nothing,
                                neighbors=nothing,
                                cutoff=nothing,
                                use_sunny_radius::Bool=true,
                                shell=:nn,
                                tol=0.2,
                                supercell=(1,1,1))

    # 1. explicit supercell
    sc, site_map = make_supercell(cryst, spec.dims)
    realized = sc
    sublog = NamedTuple[]
    vaclog = NamedTuple[]

    # 2. substitutions
    for (_, rule) in spec.substitutions
        realized, subs = dope(realized, rule; rng=rng)
        append!(sublog, subs)
    end

    # 3. vacancies
    for (_, rule) in spec.vacancies
        realized, vacs = vacate(realized, rule; rng=rng)
        append!(vaclog, vacs)
    end

    # 4. distortions
    if spec.displacement_rule !== nothing
        realized = distort(realized, spec.displacement_rule; rng=rng)
    end

    rspec = _extract_crystal(realized)

    # 5. build Phunny model
    model = Phunny.build_model(realized;
        mass=mass,
        isotopes_by_site=isotopes_by_site,
        isotopes_by_species=isotopes_by_species,
        neighbors_sunny=neighbors_sunny,
        neighbors=neighbors,
        cutoff=cutoff,
        use_sunny_radius=use_sunny_radius,
        kL=kL, kT=kT, β=β,
        shell=shell, tol=tol,
        supercell=supercell
    )

    Φ = assemble_Φ ? Phunny.assemble_force_constants!(model) : nothing

    return DisorderModel(
        parent_crystal = cryst,
        realized_crystal = realized,
        spec = spec,
        dims = spec.dims,
        fracpos = rspec.fracpos,
        species = rspec.species,
        site_map = site_map,
        substitutions = sublog,
        vacancies = vaclog,
        seed = spec.seed,
        model = model,
        Φ = Φ,
        metadata = Dict(
            :nsites => length(rspec.species),
            :nsubstitutions => length(sublog),
            :nvacancies => length(vaclog),
        )
    )
end


#---------------------------------------------------------------# 
#                       4. Ensemble Sampler                     #
#---------------------------------------------------------------#

"""
    sample_disordered_model(cryst, spec, nsamples; observable=nothing, kwargs...)

Generate `nsamples` realizations and optionally evaluate an observable
for each one.

`observable` should be a function of the form:
    observable(dm::DisorderModel) -> result

Returns `DisorderEnsemble`.
"""
function sample_disordered_model(cryst, spec::DisorderSpec, nsamples::Int;
                                 observable=nothing,
                                 rng = spec.seed === nothing ? Random.default_rng() : MersenneTwister(spec.seed),
                                 kwargs...)

    samples = DisorderModel[]
    obs = Any[]

    # use deterministic child seeds if parent seed exists
    seeds = spec.seed === nothing ? fill(nothing, nsamples) : rand(rng, 1:10^9, nsamples)

    for n in 1:nsamples
        localspec = DisorderSpec(
            dims = spec.dims,
            substitutions = spec.substitutions,
            vacancies = spec.vacancies,
            displacement_rule = spec.displacement_rule,
            seed = seeds[n],
            label = spec.label,
            metadata = copy(spec.metadata),
        )

        dm = build_disordered_model(cryst, localspec; kwargs...)
        push!(samples, dm)

        if observable !== nothing
            push!(obs, observable(dm))
        end
    end

    ens = DisorderEnsemble(
        parent_crystal = cryst,
        spec = spec,
        nsamples = nsamples,
        samples = samples,
        observables = obs,
        metadata = Dict(
            :has_observables => observable !== nothing,
        )
    )

    # simple numeric summary if possible
    if !isempty(obs)
        try
            A = stack(obs)
            ens.mean_observable = mean(A; dims=ndims(A))
            ens.var_observable  = var(A; dims=ndims(A))
        catch
            # leave summaries as nothing if object type is not stackable
        end
    end

    return ens
end
