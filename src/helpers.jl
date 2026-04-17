# --------------------------------------------#
# q-space conversion ~ (r.l.u.) --> cartesian #
# --------------------------------------------#
"""
    q_cartesian(cryst, q; basis=:cart, cell=:primitive)

Convert a wavevector to Cartesian (Å⁻¹). `q` may be:
- `basis=:cart`    → already Cartesian (Å⁻¹)
- `basis=:rlu`     → reduced lattice units (h,k,l). `cell` chooses `:primitive` or `:conventional`.
Uses `Sunny.prim_recipvecs(cryst)` if Sunny is loaded and `cell=:primitive`, otherwise falls back to
`2π * inv(L)'` where `L` is the corresponding direct lattice.
"""
function q_cartesian(cryst, q::SVector{3,Float64}; basis::Symbol=:cart, cell::Symbol=:primitive)
    if basis === :cart
        return q
    elseif basis === :rlu
        if cell === :primitive && isdefined(Main, :Sunny)
            S = _sunny_mod()
            G = S.prim_recipvecs(cryst)             # 3x3 columns b1,b2,b3 (Å⁻¹)
            return SVector{3,Float64}(G * q)
        else
            # conventional cell
            L = SMatrix{3,3,Float64,9}(Matrix(getproperty(cryst, :latvecs)))
            G = 2π * inv(L)'                         # columns b1,b2,b3
            return SVector{3,Float64}(G * q)
        end
    else
        throw(ArgumentError("basis must be :cart or :rlu"))
    end
end
#
#
#----------------------------#
# Neighbor/bond construction #
#----------------------------#
"""
    neighbor_bonds_cutoff(lattice, fracpos; cutoff, kL, kT, supercell=(1,1,1))

Construct bonds using a radial cutoff by scanning a small supercell (portable path).
`kL` and `kT` can be numbers or functions `(i,j,rij)->stiffness` for non-uniform springs.
Returns `Vector{Bond}`.
"""
function neighbor_bonds_cutoff(lattice::SMatrix{3,3,Float64,9},
                               fracpos::Vector{SVector{3,Float64}};
                               cutoff::Real,
                               kL::Union{Real,Function}=1.0,
                               kT::Union{Real,Function}=1.0,
                               supercell::NTuple{3,Int}=(1,1,1))
    S = SVector{3,Int}(supercell)
    N = length(fracpos)
    bonds = Bond{Float64}[]
    cutoff2 = cutoff^2

    for i in 1:N
        ri = frac_to_cart(lattice, fracpos[i])
        for j in 1:N
            for Rx in -(S[1]÷2):(S[1]÷2), Ry in -(S[2]÷2):(S[2]÷2), Rz in -(S[3]÷2):(S[3]÷2)
                R = SVector{3,Int}(Rx,Ry,Rz)
                rj = frac_to_cart(lattice, fracpos[j] .+ SVector{3,Float64}(R))
                rij = rj - ri
                d2 = dot(rij, rij)
                (d2 == 0.0 || d2 > cutoff2) && continue
                kL_ = kL isa Function ? kL(i,j,rij) : float(kL)
                kT_ = kT isa Function ? kT(i,j,rij) : float(kT)
                push!(bonds, Bond(i, j, R, rij, kL_, kT_))
            end
        end
    end

    # Deduplicate (canonical ordering)
    keep = trues(length(bonds))
    for (idx,b) in pairs(bonds)
        if b.i == b.j && b.R == SVector{3,Int}(0,0,0)
            keep[idx] = false; continue
        end
        if b.i > b.j || (b.i == b.j && any(b.R .< 0))
            keep[idx] = false
        end
    end
    return [bonds[i] for i in eachindex(bonds) if keep[i]]
end
"""
    neighbor_bonds_from_sunny(cryst, bonds; kL=1.0, kT=1.0)

Fast-path bond construction using Sunny's bonds and geometry.
`bonds` may be a `Vector` of:
- Sunny.Bond         (fields: i,j,n)
- Tuples (i,j,n::SVector{3,Int})
- Tuples (i,j,Δ::SVector{3,Float64})  # direct Cartesian Δ

`kL`/`kT` can be numbers or `(i,j,r0)->value`. Returns `Vector{Bond}`.
"""
function neighbor_bonds_from_sunny(cryst, bonds;
                                   kL::Union{Real,Function}=1.0,
                                   kT::Union{Real,Function}=1.0)
    L = SMatrix{3,3,Float64,9}(Matrix(getproperty(cryst, :latvecs)))
    out = Bond{Float64}[]
    S = isdefined(Main, :Sunny) ? _sunny_mod() : nothing

    for nb in bonds
        if S !== nothing && (hasproperty(nb, :i) && hasproperty(nb, :j) && hasproperty(nb, :n))
            i = getproperty(nb, :i); j = getproperty(nb, :j); n = getproperty(nb, :n)
            # Use Sunny's precise geometry
            bpos = S.BondPos(cryst, nb)
            r0   = SVector{3,Float64}(S.global_displacement(cryst, bpos))
            R    = SVector{3,Int}(n)
        elseif nb isa Tuple{Int,Int,SVector{3,Int}}
            i,j,R = nb
            # Need fractional positions for j + R - i; fetch from cryst if available
            fpos = getproperty(cryst, :positions)
            r0   = SVector{3,Float64}(L * (SVector{3,Float64}(fpos[j]) .+ SVector{3,Float64}(R) .- SVector{3,Float64}(fpos[i])))
        elseif nb isa Tuple{Int,Int,SVector{3,Float64}}
            i,j,Δ = nb
            R     = SVector{3,Int}(0,0,0)
            r0    = Δ
        else
            throw(ArgumentError("Unsupported neighbor element of type $(typeof(nb))"))
        end

        d2 = dot(r0,r0)
        (d2 == 0.0) && continue
        kL_ = kL isa Function ? kL(i,j,r0) : float(kL)
        kT_ = kT isa Function ? kT(i,j,r0) : float(kT)
        push!(out, Bond(i,j,R,r0,kL_,kT_))
    end

    # Deduplicate (canonical ordering)
    keep = trues(length(out))
    for (idx,b) in pairs(out)
        if b.i == b.j && b.R == SVector{3,Int}(0,0,0)
            keep[idx] = false; continue
        end
        if b.i > b.j || (b.i == b.j && any(b.R .< 0))
            keep[idx] = false
        end
    end
    return [out[i] for i in eachindex(out) if keep[i]]
end

"""
    neighbor_bonds_radius(cryst; cutoff, kL, kT)

Sunny-powered radius query: builds bonds using `Sunny.all_bonds_for_atom`.
Requires `using Sunny`. Returns `Vector{Bond}`.
"""
function neighbor_bonds_radius(cryst; cutoff::Real, kL::Union{Real,Function}=1.0, kT::Union{Real,Function}=1.0)
    S = _sunny_mod()
    N = length(getproperty(cryst, :positions))
    bonds_sunny = S.Bond[]
    for i in 1:N
        append!(bonds_sunny, S.all_bonds_for_atom(cryst, i, cutoff))
    end
    return neighbor_bonds_from_sunny(cryst, bonds_sunny; kL=kL, kT=kT)
end

#----------------------------------------------------------#
# Three-body Bond Angle Discovery from Canonical Bond List #
#----------------------------------------------------------#
function find_bond_angles(bonds::Vector{Bond{T}}; β::Real=zero(T), shell::Symbol=:nn, tol::Real=0.20) where T
    
    #Expand canonical bonds into directed neighbor list
    neighbor = Dict{Int, Vector{Phunny.Bond{Float64}}}()
    for b in bonds
        push!(get!(Vector{Phunny.Bond{Float64}}, neighbor, b.i), b)
        if !(b.i == b.j && b.R == SVector{3,Int}(0,0,0))
            push!(get!(Vector{Phunny.Bond{Float64}}, neighbor, b.j), Phunny.Bond{Float64}(b.j,b.i,-b.R,-b.r0,b.kL,b.kT))
        end
    end

    #For each vertex (j), select eligible neighbors and enumerate angle pairs
    angles = Angle{Float64}[]; βval = Float64(β)
    for (j, nb) in neighbor
        length(nb) < 2 && continue
        
        #Shell selection
        dists = [norm(n.r0) for n in nb]
        rmin = minimum(dists)
        sel = if shell === :nn
            [idx for (idx, d) in enumerate(dists) if d ≤ (1.0 + tol)*rmin]
        elseif shell === :all
            eachindex(nb)
        else
            error("Unknown shell criteria: $shell. Use :nn or :all.")
        end
        length(sel) < 2 && continue

        #Unique pairs with canonical ordering (atomic index, lattice offset)
        for a in 1:length(sel)-1, b in a+1:length(sel)
            ni = nb[sel[a]]; nk = nb[sel[b]]
            i, Rji = ni.j, ni.R
            k, Rjk = nk.j, nk.R

            if (i > k) || (i == k && Tuple(Rji) > Tuple(Rjk))
                i, k = k, i
                Rji, Rjk = Rjk, Rji
            end

            push!(angles, Angle(i,j,k,Rji,Rjk,βval))
        end
    end
    return angles
end

#
#
#------------------------------------------------------------------------------------------#
# Automatic look-up for ionic mass (optional: isotope mass) and coherent scattering length #
# -----------------------------------------------------------------------------------------#
#
# Parse species labels, possibly with isotope mass number (e.g., "Si-29", "O18"), also "D","T"
function _parse_species_label(x)::Tuple{Symbol,Union{Int,Nothing}}
    s = String(x); s = strip(s)
    if s == "D"; return (:H, 2) end
    if s == "T"; return (:H, 3) end
    m = match(r"^([A-Za-z]{1,2})(?:-|\s*)?(\d+)?$", s)
    if m === nothing
        return (Symbol(s), nothing)
    else
        Z = Symbol(m.captures[1])
        A = m.captures[2] === nothing ? nothing : parse(Int, m.captures[2])
        return (Z, A)
    end
end

function mass_lookup(species::AbstractVector; iso_by_site::Dict{Int,Int}=Dict{Int,Int}(), iso_by_species::Dict{Symbol,Int}=Dict{Symbol,Int}())
    n = length(species)
    out = Vector{Float64}(undef, n)
    for s in 1:n
        Zs, A_tag = _parse_species_label(species[s])
        A = haskey(iso_by_site, s) ? iso_by_site[s] :
            (A_tag !== nothing ? A_tag :
             (haskey(iso_by_species, Zs) ? iso_by_species[Zs] : nothing))
        if A !== nothing
            key = (Zs, A)
            haskey(MASS_ISO_U, key) || error("No isotopic mass for $(Zs)-$(A)")
            out[s] = MASS_ISO_U[key]
        else
            haskey(MASS_U, Zs) || error("No natural-abundance mass for element $(Zs)")
            out[s] = MASS_U[Zs]
        end
    end
    return out
end

function bcoh_lookup(species::AbstractVector; iso_by_site::Dict{Int,Int}=Dict{Int,Int}(), iso_by_species::Dict{Symbol,Int}=Dict{Symbol,Int}())
    n = length(species)
    out = Vector{Float64}(undef, n)
    @inbounds for s in 1:n
        Zs, A_tag = _parse_species_label(species[s])
        A = haskey(iso_by_site, s) ? iso_by_site[s] :
            (A_tag !== nothing ? A_tag :
             (haskey(iso_by_species, Zs) ? iso_by_species[Zs] : nothing))
        if A !== nothing
            key = (Zs, A)
            if haskey(BCOH_ISO_FM, key)
                out[s] = BCOH_ISO_FM[key]
            else
                haskey(BCOH_FM, Zs) || error("No coherent length for element $(Zs)")
                out[s] = BCOH_FM[Zs]
            end
        else
            if haskey(BCOH_FM, Zs)
                out[s] = BCOH_FM[Zs]
            else
                # Allow symbols like :O16 in the table
                if haskey(BCOH_FM, Symbol(string(Zs)))
                    out[s] = BCOH_FM[Symbol(string(Zs))]
                else
                    error("No coherent length for element $(Zs)")
                end
            end
        end
    end
    return out
end


#----------------------------#
# Model construction helpers #
#----------------------------#

# Extract Sunny-like information
@inline function _to_phunny_spec(crystal)
    if hasproperty(crystal, :latvecs) && hasproperty(crystal, :positions) && hasproperty(crystal, :types)
        L  = SMatrix{3,3,Float64,9}(Matrix(getproperty(crystal, :latvecs)))
        fp = [SVector{3,Float64}(Tuple(p)) for p in getproperty(crystal, :positions)]
        sp = Symbol.(getproperty(crystal, :types))
        return (lattice=L, positions=fp, species=sp)
    else
        L  = SMatrix{3,3,Float64,9}(getproperty(crystal, :lattice))
        fp = [SVector{3,Float64}(p) for p in getproperty(crystal, :positions)]
        sp = Symbol.(getproperty(crystal, :species))
        return (lattice=L, positions=fp, species=sp)
    end
end

"""
    build_model(crystal; mass, neighbors_sunny=nothing, neighbors=nothing, cutoff=nothing, use_sunny_radius=true, kL=1.0, kT=0.0, β=0.0, shell=:nn, tol=0.2, supercell=(1,1,1))

Create a `Model` from a Sunny-like `crystal` object. Minimal interface expected:
- Either Sunny-style fields: `crystal.latvecs` (3×3), `crystal.positions` (fractional), `crystal.types` (Strings)
- or generic fields: `crystal.lattice` (3×3), `crystal.positions` (fractional), `crystal.species` (Symbols/Strings)

Bond assembly precedence:
1. If `neighbors_sunny` (e.g. `Vector{Sunny.Bond}`) is provided → fast path via Sunny geometry.
2. Else if `neighbors` (list of tuples as documented in `neighbor_bonds_from_sunny`) is provided → fast tuple path.
3. Else if `use_sunny_radius && cutoff!=nothing && Sunny loaded` → build via `Sunny.all_bonds_for_atom`.
4. Else → fallback to portable cutoff scan `neighbor_bonds_cutoff`.
"""
function build_model(crystal; mass=:lookup, isotopes_by_site=nothing, isotopes_by_species=nothing,
                     neighbors_sunny=nothing, neighbors=nothing, cutoff=nothing, use_sunny_radius::Bool=true, 
                     kL=0.0, kT=0.0, β=0.0, shell=:nn, tol=0.2, supercell=(1,1,1))
    spec = _to_phunny_spec(crystal)
    L, fpos, species = spec.lattice, spec.positions, spec.species

    massvec = if mass === :lookup
        mass_lookup(species;
            iso_by_site = isotopes_by_site === nothing ? Dict{Int,Int}() : isotopes_by_site,
            iso_by_species = isotopes_by_species === nothing ? Dict{Symbol,Int}() : isotopes_by_species)
    elseif mass isa Dict
        [haskey(mass, s) ? Float64(mass[s]) :
         haskey(mass, String(s)) ? Float64(mass[String(s)]) :
         error("No mass for species $s in provided mass mapping") for s in species]
    else
        collect(Float64.(mass))
    end

    
    bonds = nothing
    if neighbors_sunny !== nothing
        bonds = neighbor_bonds_from_sunny(crystal, neighbors_sunny; kL=kL, kT=kT)
    elseif neighbors !== nothing
        bonds = neighbor_bonds_from_sunny((; latvecs=Matrix(L), positions=fpos), neighbors; kL=kL, kT=kT)
    elseif use_sunny_radius && cutoff !== nothing && isdefined(Main, :Sunny)
        bonds = neighbor_bonds_radius(crystal; cutoff=cutoff, kL=kL, kT=kT)
    else
        cutoff === nothing && error("Provide `cutoff` if `bonds` not given (or pass `neighbors_sunny`/`neighbors`).")
        bonds = neighbor_bonds_cutoff(L, fpos; cutoff=cutoff, kL=kL, kT=kT, supercell=supercell)
    end

    angles = find_bond_angles(bonds; β=β, shell=shell, tol=tol)

    return Model(L, fpos, species, massvec, bonds, angles, length(fpos))
end


function assign_force_constants!(model::Model, atom::Dict{Int,Symbol},
				 bonds::Dict{Tuple{Symbol,Symbol},Tuple{Float64,Float64}};
                                 angles::Dict{Tuple{Symbol,Symbol,Symbol},Float64}=Dict{Tuple{Symbol,Symbol,Symbol},Float64}())
    #2-Body Bond Stretching
    @inline canonical!(s1,s2) = s1 <= s2 ? (atom[s1], atom[s2]) : (atom[s2], atom[s1])
    for e in eachindex(model.bonds)
	    b = model.bonds[e]
	    pair = canonical!(b.i, b.j)
	    if haskey(bonds, pair)
		    kL, kT = bonds[pair]
		    model.bonds[e] = Phunny.Bond(b.i, b.j, b.R, b.r0, kL, kT)
	    end
    end

    #3-Body Bond Bending
    @inline function canonical_triplet(s1,s2,s3)
        a,c = atom[s1], atom[s3]
        a <= c ? (a, atom[s2], c) : (c, atom[s2], a)
    end
    for e in eachindex(model.angles)
        ang = model.angles[e]
        triplet = canonical_triplet(ang.i, ang.j, ang.k)
        if haskey(angles, triplet)
            model.angles[e] = Phunny.Angle(ang.i, ang.j, ang.k, ang.Rji, ang.Rjk, angles[triplet])
        end
    end
    return model
end

#----------------------------------#
# Acoustic Phonon Velocity Helpers #
#----------------------------------#

# Sound speeds along given unit directions (Cartesian Å⁻¹) via finite differences near Γ
# Returns velocities for three acoustic branches (m/s) per direction.
function sound_speeds(model::Model, Φ, dirs::Vector{SVector{3,Float64}}; cryst=nothing, dq=1e-3)
    # Unit conversion
    Åinv_to_minv = 1e10

    speeds = Vector{NTuple{3,Float64}}(undef, length(dirs))
    for (k, nhat) in pairs(dirs)
        q1 = dq * nhat                    # Å⁻¹
        E, _ = phonons(model, Φ, q1; q_basis=:cart, cryst=cryst)  # meV
        E = sort(E)
        # Approximate slope dE/dq for the three acoustic branches
        dEdq = E[1:3] ./ dq               # meV·Å
        # Convert to m/s: v = (dE/dq)/ħ  with E in J, q in m⁻¹ → multiply by (meV→J) and (Å→m)
        v = ntuple(p-> (dEdq[p] * meV_to_J / ℏ) / Åinv_to_minv, 3)
        speeds[k] = v
    end
    return speeds
end

#=
    _max_acoustic_speed(model, Φ; cryst=nothing, q_cell=:primitive, dq=1e-3)

Estimates the maximum acoustic group velocity over the three primitive reciprocal
lattice directions, returning the result in meV·Å.

Samples E(dq·n̂)/dq at a small but finite wavevector along each primitive reciprocal
direction — the direct finite-difference definition of the acoustic group velocity at Γ.
This is more numerically stable near Γ than the gradient method, which divides by ω.

Used internally to construct the q-dependent mode threshold in `onephonon_dsf`:
    eps_meV_q = eps_scale * v_max * |q_cart|

Falls back to the conventional reciprocal lattice (2π·inv(L)') if Sunny is not loaded
or no crystal is provided.

=#
function _max_acoustic_speed(model::Model,
                              Φ::Dict{Tuple{Int,Int,SVector{3,Int}},SMatrix{3,3,Float64,9}};
                              cryst=nothing, q_cell::Symbol=:primitive,
                              dq::Float64=1e-3)::Float64
    # Primitive reciprocal directions — correct for any crystal symmetry
    G = if cryst !== nothing && isdefined(Main, :Sunny) && q_cell === :primitive
        getfield(Main, :Sunny).prim_recipvecs(cryst)
    else
        2π * inv(model.lattice)'   # fallback: conventional reciprocal from direct lattice
    end

    vmax_meVA = 0.0
    for col in eachcol(G)
        nhat   = SVector{3,Float64}(col / norm(col))
        # sound_speeds returns NTuple{3,Float64} of speeds in m/s per direction
        speeds = sound_speeds(model, Φ, [nhat]; cryst=cryst, dq=dq)[1]
        # Convert m/s → meV·Å: v[meV·Å] = v[m/s] * ℏ[J·s] / (meV_to_J[J] * 1e-10[m/Å])
        for vms in speeds
            vmax_meVA = max(vmax_meVA, vms * ℏ / (meV_to_J * 1e-10))
        end
    end
    return max(vmax_meVA, 1e-2)   # 1e-2 meV·Å floor for purely optical or pathological models
end



#------------------------------------------------------#
# atomic_index() ~ maps index to atomic label	       #
#------------------------------------------------------#
"""
    atomic_index(labels::Vector{String})

Takes in atomic labels (e.g. ["Cu", "O"]) and outputs a dictionary 
which maps atomic index to atomic label (e.g. [1 => :Cu, 2 => :O])
"""
@inline atomic_index(labels) = Dict(i => Symbol(labels[i]) for i in eachindex(labels))


#------------------------------------------------------#
# make_physical!() ~ returns mass-weighted eigenvector #
#------------------------------------------------------#
@inline make_physical!(V,v,m,s) = SVector{3}(V[3(s-1)+1,v], V[3(s-1)+2,v], V[3(s-1)+3,v])/sqrt(m[s])

#-----------------------------------------------------------#
# longitudinal_weight!() ~ computes longitudinal phonon     #
# projection using mass-weighted eigenvectors & unit bond   #
# direction (nhat).					    #
#-----------------------------------------------------------#
@inline longitudinal_weight!(V,v,m,s,nhat) = abs2(dot(nhat, make_physical!(V,v,m,s)))

#------------------#
# S(q,w) Reshaping #
#------------------#
const AX = (; h=1, k=2, ℓ=3, l=3, ω=4, w=4)
"""
    collapse(A::AbstractArray; over::Symbol, op::Function)

Collapse (reduce) a 4D array `A` over one or more axes with specified operation (default operation is `sum`).
Axes may be symbols (:h,:k,:ℓ,:ω) [Julia Syntax] or integer indices [Python Syntax].

Reduces highest axes first to keep indices stable. Returns collapsed (reduced) array.
"""
function collapse(A; over=:ω, op=sum)
	axes = over isa Tuple ? over : (over,)
	idxs = sort!(map(x -> x isa Symbol ? AX[x] : Int(x), collect(axes)))
	for ax in Iterators.reverse(idxs)       # reduce highest axis first
		A = dropdims(op(A; dims=ax); dims=ax)
	end
	return A
end


#---------------#
# Utility grids #
#---------------#

ω_grid(ωmin, ωmax, n) = range(ωmin, ωmax; length=n)











