# Force-constant models, symmetry and magnon–phonon coupling

## Linear force-constant basis

Bond stretching, transverse and bond-bending stiffnesses enter the harmonic potential linearly, and so do the mass
weighting and the Fourier transform. For parameters ``k_p``

```math
\Phi = \sum_p k_p \Phi_p \quad\Longrightarrow\quad D(\mathbf q) = \sum_p k_p D_p(\mathbf q),
\qquad \frac{\partial D}{\partial k_p} = D_p .
```

[`ForceConstantBasis`](@ref) groups bonds by label pairs or by a classifier `bond -> Symbol` (which can separate
bonds by length) and angles by `angle -> Symbol`. [`dynamical_matrices`](@ref) caches the ``D_p(\mathbf q)``
once per wave vector; every further parameter set costs one weighted sum and one diagonalization.

## Fitting and posterior sampling

[`fit_force_constants`](@ref) fits positive parameters to measured values, ``\chi^2 = \sum_i (o_i(k) - d_i)^2/\sigma_i^2``,
in ``\theta = \log k`` inside a box, by multi-start Levenberg–Marquardt. The observables are any function of ``k``,
typically zone-centre frequencies from [`gamma_irreps`](@ref) and the elastic Debye temperature from
[`debye_temperature`](@ref),

```math
\Theta_D = \frac{\hbar}{k_B} v_m (6\pi^2 n)^{1/3}, \qquad v_m^{-3} = \frac13 \Big\langle \sum_{s=1}^{3} v_s^{-3} \Big\rangle_{\hat{\mathbf q}} .
```

[`sample_force_constants`](@ref) draws from the posterior ``\propto e^{-\chi^2/2}`` (uniform prior in the box) by
random-walk Metropolis, so that poorly constrained combinations of force constants are carried forward rather
than fixed at a point estimate.

## Zone-centre symmetry

[`displacement_representation`](@ref) gives the matrix of a space-group operation ``\{R|\mathbf t\}`` on the
``3N`` displacements at Γ. [`gamma_irreps`](@ref) groups degenerate modes, computes their characters under a list
of operations, and matches them to a user-supplied character table. Mass-weighted eigenvectors give the same
characters, since operations exchange only like atoms.

## Magnon–phonon hybridization

**Phases.** Phonon eigenvectors carry atom-position phases ``e^{i\mathbf q\cdot\mathbf r}``; the magnon BdG matrix
must use the same convention. The displacement of atom ``\kappa`` per unit of ``(c_{\mathbf q} + c^\dagger_{-\mathbf q})`` is
``x_0(M_\kappa, E)\, e_\kappa(\mathbf q)`` with ``x_0 = \sqrt{\hbar^2/2ME}`` ([`mode_amplitude`](@ref)).

**Local strain.** For a magnetic atom with neighbours at bond vectors ``\Delta`` (``3 \times n``), the
least-squares displacement gradient of the relative displacements ``U`` is

```math
G = U P, \qquad P = \Delta^{\mathsf T} (\Delta \Delta^{\mathsf T})^{-1},
```

exact for affine displacement fields ([`LocalEnvironment`](@ref), [`displacement_gradient`](@ref)). Its symmetric
and antisymmetric parts are the local strain ``\varepsilon`` and rotation ``\omega``, rotated into the local frame
of the moment (``\hat z'`` along the moment) of each [`MagneticSite`](@ref).

**Coupling.** A coupling is any function ``(\varepsilon', \omega', S) \mapsto (h_x, h_y)``, the fields conjugate to
``S^{x'}, S^{y'}``. To linear order ``S^{x'} = \sqrt{S/2}(a + a^\dagger)``, ``S^{y'} = -i\sqrt{S/2}(a - a^\dagger)``.
[`single_ion_magnetoelastic`](@ref) implements the uniaxial single ion,

```math
H_c = \sum_{\alpha = x', y'} \left(B\, \varepsilon_{\alpha z'} - D\, \omega_{\alpha z'}\right) O_{\alpha z'},
\qquad O_{\alpha z'} \to (2S - 1) S_\alpha ,
```

whose rotational term is required by rotational invariance (Melcher 1972).

**Spectrum.** [`magnon_phonon_bdg`](@ref) assembles ``M = \begin{pmatrix} A(\mathbf q) & B(\mathbf q) \\ B(\mathbf q)^\dagger & A(-\mathbf q)^* \end{pmatrix}``
in the basis ``(a_{\mathbf q}, c_{\mathbf q}, a^\dagger_{-\mathbf q}, c^\dagger_{-\mathbf q})``, with ``e(-\mathbf q) = e(\mathbf q)^*``.
[`colpa`](@ref) diagonalizes it (Colpa 1978): ``M = K^\dagger K``, then the Hermitian problem of ``K \sigma_3 K^\dagger``,
giving a para-unitary ``T``. A failed Cholesky factorization means the assumed ordered state is unstable.
[`hybrid_modes`](@ref) returns the magnon fraction of each mode, and [`first_order_coupling`](@ref) the matrix
elements between uncoupled magnons and phonons, so that a crossing detuned by ``\delta`` splits by
``\sqrt{\delta^2 + 4|V|^2}``.

## References

- J. H. P. Colpa, *Physica A* **93**, 327 (1978).
- R. L. Melcher, *Phys. Rev. Lett.* **28**, 165 (1972).
