using LinearAlgebra, StaticArrays, Sunny, GLMakie, Statistics, BenchmarkTools

include(joinpath(@__DIR__, "../..", "src/Phunny.jl"))
using .Phunny

#---------------------------------------------------------------#
#References: 							#
#	     - Jorgensen et al., Phys. Rev. B 38, 11337 (1988) 	#
#	     - Radaelli et al., Phys. Rev. B 49, 4163 (1994) 	#
#---------------------------------------------------------------#



################
# Define Model #
################


#Lattice Geometry : Low Temperature, Space Group: Bmab
a,b,c = 5.35/sqrt(2), 5.41/sqrt(2), 13.2 
L = lattice_vectors(a, b, c, 90, 90, 90)

#Small in-plane tilt applied to planar oxygens
δ = 0.0

#Fractional Positions + Atomic Labels
fpos = [@SVector[0.0,        0.0,      0.0  ], # Cu
	@SVector[0.5 + δ,    0.0,      0.0  ], # O ~ planar (x)
	@SVector[0.0,      0.5 + δ,    0.0  ], # O ~ planar (y)
	@SVector[0.0,        0.0,      0.183], # O ~ apical (+z)
	@SVector[0.0,        0.0,      0.817], # O ~ apical (-z)
	@SVector[0.0,        0.0,      0.361], # La (+z),
	@SVector[0.0,        0.0,      0.639]] # La (-z)
types = ["Cu", "O", "O", "O", "O", "La", "La"]

#Crystal + Model
cryst = Crystal(L, fpos; types) 
cutoff = 0.75*minimum((a,b,c)); 
model = build_model(cryst; cutoff, β=0.2, shell=:all, tol=0.15)

#Manually Defined Bonds & Force Constants ~ (:Atom₁, :Atom₂) => (kL, kT)
bond_dict = Dict( (:Cu, :O)  => (10.0, 0.0),
		  (:O,  :O)  => (2.00, 0.0),
		  (:O,  :La) => (5.00, 0.0))

#Index atoms
atoms = atomic_index(types)

#Mass Sanity Check
@info "Mass (Min, Max) : $(extrema(model.mass))" 


############
# Analysis #
############
#Calibration
for (kL_CuO, kT_CuO) in ((1,0), (4,0), (10,0), (33,0), (50,0), (80, 0), (190,0), (210,0), (220, 0))
    bond_dict[:Cu,:O] = (kL_CuO, kT_CuO)
    
    assign_force_constants!(model, atoms, bond_dict)
    FCMs = assemble_force_constants!(model)

    EΓ, _ = phonons(model, FCMs, @SVector[0.0, 0.0, 0.0]; 
                    q_basis=:rlu, q_cell=:primitive, cryst=cryst)
    if any( ω -> 60.0 ≤ ω ≤ 95.0, EΓ)
        @info "Cu-O (kL=$(kL_CuO)) ⟶ Γ optics ≈  $(round.(EΓ[4:end];digits=2))"
    end
end
bond_dict[:Cu,:O] = (10.0, 0.0)
#kL : Controls Acoustic Scaling
#kT : Controls Optical Scaling

#################
# FCM & Phonons #
#################

#Assign force constants using bond dictionary
assign_force_constants!(model, atoms, bond_dict)

#Assemble Force Constants Matrices
ϕ = assemble_force_constants!(model)

#Minimum Image
@inline minimg(x) = x .- round.(x)

#Convert Fractional -> Cartesian (Real-Space)
@inline r_cartesian(Lat, r1, r2) = norm(Lat*minimg(r2 - r1))

planar_x = r_cartesian(L, fpos[1], fpos[2])
planar_y = r_cartesian(L, fpos[1], fpos[3])
apical = r_cartesian(L, fpos[1], fpos[4])

@assert 1.85 - δ ≤ planar_x ≤ 1.95 + δ && 1.85 - δ ≤ planar_y ≤ 1.95 + δ "Planar oxygen bond distance does NOT match the experimental values!"
@assert 2.30 ≤ apical ≤ 2.45 "Apical oxygen bond distance does NOT match the experimental values!"

#DSF-based check: stable phonons & sensible optical energies at Γ-point
eigenpairs = phonons(model, ϕ, @SVector[0.0, 0.0, 0.0]; 
		     q_basis=:rlu, q_cell=:primitive, cryst=cryst)

@assert any(ω -> 60.0 ≤ ω ≤ 95.0, eigenpairs[1]) "Expected oxygen-dominant optical modes in ~[60, 95] meV window for La2CuO4! $(eigenpairs[1])"
#@show size(eigenpairs[2])
#display(heatmap(abs.(eigenpairs[2]).^2))

############
# Plotting #
############


#Plot Phonon DoS
#display(hist(eigenpairs[1], bins=20, color="royalblue", label="Phonon DoS"))


#Plots S(q,ω) = ∑ₙ S(qₙ, ω) : S(qₙ,ω) = S[qₙ,:]
function plot_dsf_line!(cryst, model, Φ; 
                        q₀=@SVector[-2.0, 0.0, 0.0], 
                        q₁=@SVector[2.0, 0.0, 0.0], 
                        nq=101, ωmax=110.0, nω=801, 
                        η=1.0, q_cell=:primitive)
    qs = [SVector{3,Float64}((1-t).*q₀ .+ t.*q₁) for t in range(0,1;length=nq)]
    ωs = range(0.0, ωmax; length=nω)
    σ = η/sqrt(8*log(2))

    # Precompute anisotropic Debye–Waller tensors once (Å^2)
    U₁ = U_from_phonons(model, Φ; T=0.0, cryst=cryst, qgrid=(16,16,16), q_cell=q_cell)
    U₂ = U_from_phonons(model, Φ; T=600.0, cryst=cryst, qgrid=(16,16,16), q_cell=q_cell)
    
    Sqω₁ = zeros(Float64, nq, nω)
    Sqω₂ = zeros(Float64, nq, nω)
    @inbounds for (iq, q) in enumerate(qs)
        Sω₁ = onephonon_dsf(model, Φ, q, ωs; q_basis=:rlu, q_cell=q_cell, cryst=cryst, T=0.0, _U_internal=U₁)
        Sω₂ = onephonon_dsf(model, Φ, q, ωs; q_basis=:rlu, q_cell=q_cell, cryst=cryst, T=600.0, _U_internal=U₂)
        Sqω₁[iq,:] .= Sω₁
        Sqω₂[iq,:] .= Sω₂
    end
    
    @inline lohi(z) = begin
	lo = 0.0
	m = (isfinite.(z) .& ( z .> 0.0))
	hi = any(m) ? mean(z[m]) : maximum(z)/2
	(lo, hi)
    end

    fig = Figure(size=(800,400))
    
    ax₁ = Axis(fig[1,1], xlabel="q index (X ↦ Γ ↦ X)", ylabel = "Energy (meV)", title="La₂CuO₄ | One-Phonon S(q,ω) | T = 0K") 
    hm = heatmap!(ax₁, 1:nq, ωs, Sqω₁ ; interpolate=true, colormap=:viridis, colorrange=lohi(Sqω₂))
    Colorbar(fig[1,2], hm; label="Intensity")
    
    ax₂ = Axis(fig[1,3], xlabel="q index (X ↦ Γ ↦ X)", ylabel = "Energy (meV)", title="La₂CuO₄ | One-Phonon S(q,ω) | T = 600K") 
    hm = heatmap!(ax₂, 1:nq, ωs, Sqω₂ ; interpolate=true, colormap=:viridis, colorrange=lohi(Sqω₂))
    Colorbar(fig[1,4], hm; label="Intensity") 

    screen=display(fig); wait(screen)
end
Base.@time plot_dsf_line!(cryst, model, ϕ; q_cell=:primitive, η=0.5, nq=601, nω=1001)

velocities = group_velocity(model, ϕ, @SVector[1.0, 0.0, 0.0]; cryst=cryst, nhat=nothing, dq=1e-3)

for (phonon_mode, velocity) in enumerate(velocities)
    velocity ≤ 0.1 && continue
    @show (phonon_mode, velocity)
end

summary = validate_summary(model, ϕ; cryst=cryst, qpath_rlu=nothing, T=0.0)
@show summary[:ASR_residual]
@show summary[:Rigid_translation_residual]
@show summary[:Rigid_rotation_force_net] 
@show summary[:Rigid_rotation_torque_net] 
@show summary[:Rigid_rotation_force_max]



