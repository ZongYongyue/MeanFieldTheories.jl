# Collective Excitations on Top of the Mean Field

---

## 1. The Physical Problem

A central question in many-body physics is: given a system whose ground state is known (or approximately known), what are its low-energy excitations? These excitations — phonons, magnons, excitons, magnetorotons — govern the dynamical response of the system and determine virtually all experimentally measurable spectral properties.

In many settings, the ground state is obtained within a mean-field framework such as Hartree-Fock (HF). The mean-field theory provides a set of quasiparticle bands $E^n_\mathbf{k}$ and a ground state $|G\rangle$ constructed by filling the occupied bands. The simplest excitations one can imagine are free particle-hole pairs: promote an electron from an occupied band $n_0$ to an unoccupied band $n$, costing energy $E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}$. This is the mean-field picture of excitations — a continuum of independent particle-hole pairs.

But this picture misses an essential piece of physics: **the residual interactions between particles and holes**. These interactions can bind particle-hole pairs into collective modes whose energies and dispersions are qualitatively different from the free particle-hole continuum. For example, magnons in an antiferromagnet have a dispersion that goes to zero at long wavelengths (Goldstone's theorem), while the free particle-hole continuum has a finite gap everywhere. To capture such collective excitations, one must go beyond the mean-field description of the excitation spectrum, even if the ground state itself remains at the mean-field level.

This document traces three approaches to this problem — the single-mode approximation, the particle-hole variational theory, and the Bethe-Salpeter equation — and explains how they are related.

---

## 2. Feynman's Single-Mode Approximation (1954)

### 2.1 Historical Context

In the 1940s, Landau proposed the two-fluid model for superfluid $^4\text{He}$ and postulated an energy-momentum dispersion relation for its elementary excitations: linear (phonon-like) at small momenta, with a minimum (the "roton minimum") at a finite momentum. This dispersion successfully explained the thermodynamic properties of the superfluid, but it was purely phenomenological — not derived from microscopic theory.

In 1954, Richard Feynman provided the first microscopic derivation of this dispersion curve, using a remarkably simple argument based on the variational principle. This work, published as "Atomic Theory of the Two-Fluid Model of Liquid Helium" (Physical Review 94, 262), introduced what is now called the **single-mode approximation (SMA)**.

### 2.2 The Construction

Feynman's key insight was that in a Bose liquid, due to the symmetry of the wavefunction under particle exchange, there are no low-lying single-particle excitations. The only low-energy excitations are **collective density oscillations**. He therefore proposed the following family of excited states, labeled by wavevector $\mathbf{q}$:

```math
|\mathbf{q}\rangle = \rho_\mathbf{q}\, |G\rangle, \qquad \rho_\mathbf{q} = \sum_j e^{i\mathbf{q}\cdot\mathbf{r}_j}
```

where $\rho_\mathbf{q}$ is the Fourier component of the density operator and $|G\rangle$ is the ground state. For $\mathbf{q} \neq 0$, this state is automatically orthogonal to $|G\rangle$ because a uniform liquid has no density fluctuations in its ground state.

A helpful analogy, due to Girvin: for a harmonic oscillator with ground state $\psi_0(x) \propto e^{-\alpha x^2/2}$, the first excited state is simply $\psi_1(x) = x\,\psi_0(x)$ — multiplying by $x$ flips the parity, ensuring orthogonality to $\psi_0$, and happens to give the exact first excited state. Feynman's $\rho_\mathbf{q}$ plays the same role as $x$: it is the simplest operator that creates a density modulation at wavevector $\mathbf{q}$.

Crucially, for each $\mathbf{q}$ the excited state is **completely determined with no adjustable parameters**. Feynman did not minimize over a parameter space. Instead, the variational principle serves as a **guarantee of an upper bound**: $\langle \mathbf{q}|H|\mathbf{q}\rangle / \langle \mathbf{q}|\mathbf{q}\rangle \geq \varepsilon_{\text{true}}(\mathbf{q})$, where $\varepsilon_{\text{true}}(\mathbf{q})$ is the true lowest excitation energy in the momentum-$\mathbf{q}$ sector. The tightness of this bound depends on how well a single density wave captures the actual lowest excitation.

### 2.3 The Result

The excitation energy evaluates to the celebrated formula:

```math
\varepsilon(\mathbf{q}) = \frac{f(\mathbf{q})}{S(\mathbf{q})}
```

where $f(\mathbf{q}) = \hbar^2 q^2 / 2m$ is the $f$-sum rule (energy-weighted sum rule of the dynamic structure factor), and $S(\mathbf{q}) = \langle G|\rho_{-\mathbf{q}}\rho_\mathbf{q}|G\rangle / N$ is the **static structure factor** of the ground state.

The physical content of this formula is striking: **the excitation spectrum is entirely determined by a ground-state property**. At long wavelengths, $S(q) \propto q$ for a compressible liquid, giving the linear phonon dispersion $\varepsilon(q) \propto q$. At finite $q$, the peak in $S(\mathbf{q})$ (arising from short-range correlations at the interparticle spacing) produces the roton minimum — exactly as Landau had postulated.

### 2.4 Why "Single Mode"?

The name reflects the core assumption: at each wavevector $\mathbf{q}$, there is only **one dominant collective mode** (the density wave) that carries the overwhelming majority of the spectral weight. If the dynamic structure factor $S(\mathbf{q}, \omega)$ is sharply peaked at a single frequency, then $f(\mathbf{q})/S(\mathbf{q})$ gives the position of that peak. If spectral weight is spread over a broad continuum, the SMA gives only an average — and the approximation degrades.

### 2.5 Feynman-Cohen Improvement (1956)

Feynman and Cohen later improved the theory by incorporating **backflow corrections**: the excited-state ansatz was modified to include the response of surrounding particles to the density perturbation. This improved the quantitative prediction of the roton gap, bringing it into much better agreement with neutron scattering experiments.

---

## 3. Extension to the Fractional Quantum Hall Effect: Girvin-MacDonald-Platzman (1985–86)

### 3.1 The GMP Theory

In a seminal work, Girvin, MacDonald, and Platzman (GMP) extended Feynman's SMA to the **fractional quantum Hall effect (FQHE)**, developing the **magneto-roton theory** of collective excitations (Physical Review B 33, 2481, 1986).

The strategy directly parallels Feynman's: take the Laughlin wavefunction as the known ground state $|G\rangle$, and act on it with the density operator $\bar{\rho}_\mathbf{q}$ projected into the lowest Landau level. The excitation energy is again given by the ratio of the $f$-sum rule to the static structure factor.

### 3.2 Key Results

The GMP theory made several important predictions:

- **Incompressibility implies a gap**: Unlike superfluid $^4\text{He}$ where the phonon mode is gapless, FQHE states are incompressible. Their structure factor behaves as $S(q) \propto q^4$ at long wavelengths (rather than $q$), leading to a **finite energy gap** for collective excitations at $q = 0$.
- **Magneto-roton minimum**: At finite $q$ (near the reciprocal lattice vector of the corresponding Wigner crystal), the dispersion develops a minimum — the magneto-roton minimum, analogous to the roton in helium.
- **Quantitative accuracy**: The SMA predictions agree remarkably well with exact diagonalization results for small systems.

### 3.3 The GMP Algebra

A deep insight from the GMP work is that density operators projected into a single Landau level do not commute. They satisfy a nontrivial algebraic structure now known as the **GMP algebra**, which plays a fundamental role in the modern understanding of the topology and geometry of quantum Hall states.

### 3.4 Haldane: Pseudopotentials, Model Wavefunctions, and Geometry

Haldane's contributions deepened the SMA/GMP framework in several important ways:

- **Haldane pseudopotentials**: Haldane developed a pseudopotential framework for analyzing many-body interactions in the lowest Landau level, decomposing the interaction into components labeled by the relative angular momentum of particle pairs. This provided the mathematical language in which the GMP SMA results could be most cleanly understood, and established criteria for when Laughlin-type ground states are exact.

- **Exact diagonalization benchmarks**: The early exact diagonalization studies by Haldane and Rezayi on small systems (on the sphere and torus geometries) provided the numerical benchmarks against which the GMP SMA predictions were tested — and found to be in remarkable quantitative agreement.

- **Model wavefunctions for collective modes (2012)**: In a collaboration with Yang, Hu, and Papić (Physical Review Letters 108, 256807, 2012), Haldane and coworkers went beyond the GMP approach by constructing **explicit model wavefunctions** for the magneto-roton mode on the sphere. The original GMP theory gave the excitation *energy* via the SMA but did not provide the excited-state *wavefunction* itself. By using the SMA operator projected onto the spherical geometry and analyzing its action on the Laughlin state, they obtained concrete many-body wavefunctions for the neutral collective excitations, which could be directly compared with exact diagonalization eigenstates. This work also extended the magneto-roton theory to non-Abelian states such as the Moore-Read state.

- **Geometric description and the chiral graviton mode**: Haldane proposed a geometric interpretation of the long-wavelength ($q \to 0$) limit of the GMP collective mode. In this picture, the neutral excitation at $q = 0$ reflects an oscillation of the **intrinsic metric** of the Landau orbitals — a quadrupolar deformation of the correlation hole surrounding each electron. This spin-2 excitation is now called the **chiral graviton mode (CGM)**. It provides a fundamentally new perspective on what the GMP density mode *is* at long wavelengths: not just a density fluctuation, but a fluctuation of the internal geometry of the quantum Hall fluid. The chiral graviton was experimentally observed in 2024 (Liang et al., Nature 628, 78).

### 3.5 From Continuum to Lattice: Fractional Chern Insulators

The discovery of **fractional Chern insulators (FCIs)** — lattice analogs of FQHE states that arise in topological flat bands without an external magnetic field — raised the question of whether the SMA and magneto-roton physics carry over to the lattice setting.

Repellin et al. (Physical Review B 90, 045114, 2014) systematically tested the SMA for FCIs on several lattice models, including the Haldane honeycomb model, the kagome model, and the ruby lattice model. They found that:

- The magneto-roton mode can be identified in FCIs and contains the same number of states as in the continuum FQH case, provided an appropriate mapping of momentum quantum numbers is applied.
- The SMA captures the dispersive magneto-roton branch well for "good" FCI models (those with flat bands, small ground-state splitting, and a clear entanglement gap).
- The quantitative agreement between FQH and FCI magneto-roton spectra is remarkable — for instance, the many-body gap extrapolates to nearly the same value in the thermodynamic limit.

The theoretical underpinning of this correspondence is that the algebra of band-projected density operators in a Chern band reduces to the GMP algebra in the long-wavelength and constant-Berry-curvature limit. This algebraic isomorphism is ultimately why the SMA works for FCIs.

---

## 4. The Full Particle-Hole Variational Theory

### 4.1 Beyond the Fixed Operator

The SMA fixes the operator (e.g., $\rho_\mathbf{q}$) and obtains one energy per wavevector. A natural question arises: **what if we do not fix the operator, but instead variationally optimize it?**

This leads to the particle-hole variational theory. Instead of a single predetermined density-wave state, we parametrize the excited state as a **general superposition** of all possible particle-hole pairs:

```math
|\mu, \mathbf{q}\rangle = \sum_{\mathbf{k}, n} \psi^n_\mathbf{k}(\mu, \mathbf{q})\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle
```

where $n_0$ is the occupied (hole) band, $n$ runs over unoccupied (particle) bands, $\mathbf{k}$ is the hole momentum, and $\psi^n_\mathbf{k}(\mu, \mathbf{q})$ is the **envelope function** — a genuine variational parameter to be determined.

### 4.2 Comparison with SMA

| Aspect | SMA | Particle-hole variational theory |
|:-------|:----|:---------------------------------|
| Excited-state form | $\hat{A}_\mathbf{q}\|G\rangle$ with fixed operator $\hat{A}$ | $\sum_{\mathbf{k},n} \psi^n_\mathbf{k}\, f^\dagger_{\mathbf{k}+\mathbf{q},n} f_{\mathbf{k},n_0}\|G\rangle$ with variational $\psi$ |
| Number of modes per $\mathbf{q}$ | One | Multiple ($\mu = 1, 2, 3, \ldots$) |
| Adjustable parameters | None | The full envelope function $\psi^n_\mathbf{k}$ |
| Optimization | No optimization; energy is an upper bound | Rayleigh quotient minimization $\to$ eigenvalue problem |
| Result | One energy $\varepsilon(\mathbf{q})$ per wavevector | A complete excitation spectrum $\varepsilon_\mu(\mathbf{q})$ |

The SMA can be viewed as a **special case** of the particle-hole variational theory where the envelope function is pre-determined (frozen to the form dictated by the chosen operator), with no variational freedom remaining.

### 4.3 The Eigenvalue Problem

Since the envelope function $\psi^n_\mathbf{k}$ is a true variational parameter, the excitation energy becomes a Hermitian quadratic form in $\psi$. Minimizing the Rayleigh quotient under the normalization constraint leads to a standard eigenvalue problem:

```math
\sum_{\mathbf{p}, n'} \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\, \psi^{n'}_\mathbf{p}(\mu, \mathbf{q}) = \varepsilon_\mu(\mathbf{q})\, \psi^n_\mathbf{k}(\mu, \mathbf{q})
```

The effective Hamiltonian $\mathcal{H}^{\text{eff}}$ decomposes into:

- A **diagonal part**: the free particle-hole pair energy $E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}$, coming entirely from the mean field.
- **Off-diagonal parts**: the direct and exchange particle-hole interaction kernels $\mathcal{K}^d$ and $\mathcal{K}^x$, which are residual interactions not captured by the mean field.

It is precisely the off-diagonal interaction kernels that bind particle-hole pairs into collective modes. Without them, the eigenvalues would simply be the free particle-hole energies — the mean-field excitation spectrum.

### 4.4 What This Method Captures

Diagonalizing $\mathcal{H}^{\text{eff}}(\mathbf{q})$ at each $\mathbf{q}$ yields:

- The **eigenvalues** $\varepsilon_\mu(\mathbf{q})$: the complete dispersion relation of collective excitations (multiple branches labeled by $\mu$).
- The **eigenvectors** $\psi^n_\mathbf{k}(\mu, \mathbf{q})$: the internal wavefunction of each collective mode, describing how it is composed of elementary particle-hole pairs.

This is a qualitative upgrade over SMA: instead of one energy per wavevector, we obtain the **full excitation spectrum** within the single particle-hole pair subspace.

---

## 5. Connection to the Bethe-Salpeter Equation

### 5.1 An Independent Historical Thread

The eigenvalue problem derived in §4 is, mathematically, identical to a well-known equation in quantum field theory: the **Bethe-Salpeter equation (BSE)** in the **Tamm-Dancoff approximation (TDA)**.

The BSE was introduced by Salpeter and Bethe in **1951** (Physical Review 84, 1232) — three years *before* Feynman's SMA — in the context of **relativistic quantum field theory**. Their original goal was to describe two-particle bound states (e.g., positronium) by summing ladder diagrams in the particle-particle or particle-hole channel. The BSE is essentially a Dyson equation for the two-body Green's function.

The BSE and SMA thus have **completely independent origins**:

- **BSE (1951)**: relativistic QFT, two-body bound states, diagrammatic perturbation theory.
- **SMA (1954)**: superfluid $^4\text{He}$, variational principle, physical intuition about collective density oscillations.

They converged in condensed matter physics decades later, when it was recognized that the particle-hole variational theory (the natural generalization of SMA) is mathematically the same as the BSE applied to the particle-hole channel.

### 5.2 TDA vs. Full RPA

The eigenvalue problem we derived includes only particle-hole ($ph$) excitations — the state $f^\dagger_{\mathbf{k}+\mathbf{q},n} f_{\mathbf{k},n_0}|G\rangle$ creates a particle above and a hole below the Fermi level. This corresponds to the **Tamm-Dancoff approximation (TDA)**.

The full **random phase approximation (RPA)**, or equivalently the full BSE, additionally includes hole-particle ($hp$) backward excitations ($f^\dagger_{\mathbf{k},n_0} f_{\mathbf{k}+\mathbf{q},n}|G\rangle$, which "de-excite" the ground state). This leads to a non-Hermitian eigenvalue problem with the $2 \times 2$ block structure:

```math
\begin{pmatrix} A & B \\ -B^* & -A^* \end{pmatrix} \begin{pmatrix} X \\ Y \end{pmatrix} = \varepsilon \begin{pmatrix} X \\ Y \end{pmatrix}
```

where $A = \mathcal{H}^{\text{eff}}$ is the TDA matrix and $B$ encodes ground-state correlations. TDA corresponds to setting $B = 0$.

### 5.3 Terminology Across Fields

The same mathematical structure appears under different names in different communities:

| Community | Typical name | Notes |
|:----------|:-------------|:------|
| Quantum chemistry | Tamm-Dancoff approximation (TDA) | Used in TD-DFT and CI-singles for molecular excited states |
| Semiconductor / exciton physics | Bethe-Salpeter equation (BSE) | Often built on top of GW quasiparticles |
| Quantum magnetism / spin waves | Particle-hole excitation theory | Sometimes called "generalized SMA" or just "magnon calculation" |
| Nuclear physics | Tamm-Dancoff / RPA | One of the earliest applications of these methods |

These are all the same equation. The choice of terminology is largely a matter of community convention rather than physics.

---

## 6. The Hierarchy of Approximations

The three approaches discussed form a clear hierarchy:

### Level 1: Single-Mode Approximation (SMA)

- **Operator**: fixed (e.g., density operator $\rho_\mathbf{q}$)
- **Variational freedom**: none
- **Output**: one energy per wavevector (upper bound)
- **Computational cost**: minimal (requires only the ground-state static structure factor)
- **When it works well**: when spectral weight is concentrated in a single collective mode

### Level 2: Particle-Hole Variational Theory (TDA / BSE without $hp$)

- **Operator**: optimized variationally within the single particle-hole pair subspace
- **Variational freedom**: the full envelope function $\psi^n_\mathbf{k}$
- **Output**: complete excitation spectrum $\varepsilon_\mu(\mathbf{q})$ with multiple branches
- **Computational cost**: diagonalization of $\mathcal{H}^{\text{eff}}(\mathbf{q})$ for each $\mathbf{q}$
- **When it works well**: when single particle-hole pair excitations dominate (no strong ground-state correlations)

### Level 3: Full RPA / BSE with $hp$ (ground-state correlations included)

- **Operator**: optimized over both $ph$ and $hp$ channels
- **Variational freedom**: both forward ($X$) and backward ($Y$) amplitudes
- **Output**: complete spectrum; Goldstone modes are guaranteed for spontaneously broken continuous symmetries
- **Computational cost**: diagonalization of a non-Hermitian matrix of twice the dimension
- **When it matters**: when ground-state correlations are significant (the $B$ matrix is non-negligible)

The SMA sits at one extreme (maximum simplicity, minimum variational freedom), and the full RPA/BSE sits at the other. The particle-hole variational theory (TDA) occupies the middle ground — it is the natural generalization of SMA that retains a Hermitian eigenvalue problem while providing a complete excitation spectrum within the single particle-hole subspace.

---

## 7. Applications and Modern Developments

### 7.1 Quantum Magnets and Spin Waves

In magnetic systems, the particle-hole excitation corresponds to a **magnon** (spin flip). The AKLT model (Affleck-Kennedy-Lieb-Tasaki, 1987) is a classic example where the SMA gives an analytical estimate of the Haldane gap in integer spin chains. The full particle-hole variational theory (Level 2) gives the complete magnon band structure, which can be compared directly with inelastic neutron scattering experiments.

### 7.2 Tensor Network Methods

When the ground state is represented as a matrix product state (MPS) or projected entangled pair state (PEPS), the SMA idea naturally integrates into the tensor network framework: replace one tensor in the network with a different tensor, then Fourier transform to obtain a momentum eigenstate. This has been formalized as the **MPS excitation ansatz** (Haegeman et al., 2013), which can be viewed as a tensor-network realization of the particle-hole variational theory — the replacement tensor is variationally optimized, going beyond the fixed-operator SMA.

### 7.3 Quantum Monte Carlo

The SMA formula $\varepsilon(\mathbf{q}) = f(\mathbf{q})/S(\mathbf{q})$ requires only equal-time correlation functions (the static structure factor), which are efficiently computable in quantum Monte Carlo (QMC) simulations. This makes SMA a practical route to extract excitation spectra from QMC ground-state data, which otherwise works in imaginary time and has difficulty accessing real-frequency dynamics directly.

### 7.4 Excitons and Optical Response in Solids

In semiconductor physics, the BSE (Level 2 or 3, typically built on top of GW quasiparticles rather than Hartree-Fock) is the standard method for computing **excitonic spectra** and **optical absorption**. The particle-hole interaction kernels $\mathcal{K}^d$ and $\mathcal{K}^x$ correspond to the screened and bare Coulomb interactions between the electron and hole forming the exciton. This is perhaps the most widespread application of the BSE in modern computational condensed matter physics.

### 7.5 Moiré Systems and Flat-Band Materials

In recent studies of moiré materials (e.g., twisted MoTe$_2$, twisted bilayer graphene), the BSE/particle-hole variational theory has been applied on top of Hartree-Fock ground states to compute magnon and exciton dispersions in flat-band systems with strong correlations. The interplay of topology, valley quantum numbers, and interactions makes the structure of the direct and exchange kernels particularly rich in these systems.

---

## 8. Summary

The problem of collective excitations on top of a mean-field ground state has been approached from multiple directions over the past seven decades:

- **Feynman (1954)** introduced the **single-mode approximation** for superfluid helium, using a fixed density operator to construct excitations and deriving the celebrated $\varepsilon = f/S$ formula.
- **Girvin, MacDonald, and Platzman (1985–86)** extended the SMA to the fractional quantum Hall effect, predicting the magneto-roton dispersion.
- **Salpeter and Bethe (1951)**, working independently in the context of relativistic QFT, derived the **Bethe-Salpeter equation** for two-particle bound states.
- The **particle-hole variational theory** — where the envelope function is a true variational parameter, leading to an eigenvalue problem — is the natural generalization of the SMA. It is mathematically identical to the BSE in the **Tamm-Dancoff approximation**.

These three threads, originating from superfluid physics, quantum field theory, and variational many-body theory, converge to the same mathematical structure: **diagonalizing an effective Hamiltonian in the particle-hole pair space**. The eigenvalues give the collective excitation spectrum; the eigenvectors give the internal structure of each collective mode.
