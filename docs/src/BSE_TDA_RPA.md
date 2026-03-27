# Bethe-Salpeter Equation, Tamm-Dancoff and Random Phase Approximation

---

## 1. The Bethe-Salpeter Equation

The **Bethe-Salpeter equation (BSE)** is the exact equation of motion for the two-particle Green's function $L$. It was introduced by Salpeter and Bethe in 1951 (Physical Review 84, 1232) in the context of relativistic bound states, and later became a cornerstone of many-body perturbation theory in condensed matter physics.

### 1.1 Two-Particle Green's Function

The central object is the four-point correlation function (two-particle Green's function):

```math
L(1,2;1',2') = G_2(1,2;1',2') - G(1,1')\,G(2,2')
```

where $G_2$ is the full two-particle Green's function and $G$ is the one-particle Green's function. $L$ describes the correlated propagation of two particles (or a particle-hole pair) beyond what is already captured by independent propagation.

### 1.2 The Exact BSE

The BSE is a Dyson-like equation for $L$:

```math
L = L_0 + L_0\, \Xi\, L
```

where:
- $L_0 = G\,G$ is the **non-interacting two-particle propagator** (a product of two dressed one-particle Green's functions)
- $\Xi$ is the **irreducible particle-hole interaction kernel** — the sum of all two-particle irreducible diagrams in the particle-hole channel

This equation is **exact**: no approximation has been made. The approximation enters through the choice of $G$ and $\Xi$.

### 1.3 What Makes BSE Exact (and Hard)

The exact BSE requires:

1. **Fully dressed propagators**: $G$ should be the exact interacting single-particle Green's function, including all self-energy corrections $\Sigma$ to infinite order.

2. **Irreducible kernel $\Xi$**: must include all diagrams that are two-particle irreducible in the particle-hole channel — not just the bare interaction, but also vertex corrections, screening, and higher-order processes.

3. **Frequency dependence**: both $G$ and $\Xi$ are frequency-dependent, making the BSE an integral equation in both momentum and frequency.

In practice, no one solves the exact BSE. Instead, different levels of approximation to $G$ and $\Xi$ define a hierarchy of methods.

---

## 2. Hierarchy of Approximations

### 2.1 Level 0: Hartree-Fock + Bare Interaction

Replace:
- $G \to G_{\text{HF}}$: the Hartree-Fock Green's function (single-particle energies and wavefunctions from SCF)
- $\Xi \to V_{\text{bare}}$: the bare (unscreened) Coulomb or Hubbard interaction

This is the level implemented in **MeanFieldTheories.jl**. It gives the standard **RPA** and **TDA** of condensed matter and nuclear physics.

### 2.2 Level 1: GW + Screened Interaction

Replace:
- $G \to G_{\text{GW}}$: quasiparticle Green's function with self-energy $\Sigma = iGW$ (the GW approximation)
- $\Xi \to W$: the screened Coulomb interaction (within RPA screening)

This is the standard **GW-BSE** approach widely used in computational materials science for excitonic spectra and optical absorption. The screening of the interaction and the self-energy correction to the band gap are both crucial for quantitative accuracy in real materials.

### 2.3 Level 2: Beyond GW

Include vertex corrections in both $\Sigma$ and $\Xi$:
- $\Sigma$ includes vertex corrections beyond GW (e.g., the GWΓ approximation)
- $\Xi$ includes diagrams beyond the screened interaction (e.g., second-order exchange, T-matrix contributions)

This level is largely a frontier of current research and is rarely attempted in practice.

### 2.4 Summary Table

| Level | Propagator $G$ | Kernel $\Xi$ | Method name |
|:------|:---------------|:-------------|:------------|
| Exact | Full $G$ | Full irreducible $\Xi$ | Exact BSE |
| 2 | GWΓ | $W$ + vertex corrections | Beyond-GW BSE |
| 1 | GW | Screened $W$ | GW-BSE |
| **0** | **Hartree-Fock** | **Bare $V$** | **HF-RPA / HF-TDA** |

Going down the table trades accuracy for computational simplicity. MeanFieldTheories.jl operates at Level 0, which is appropriate for model Hamiltonians (Hubbard, Heisenberg-like) where the interaction is already short-ranged and the primary interest is in qualitative collective mode structure (magnon dispersions, Goldstone modes, excitation gaps).

---

## 3. TDA and RPA: Two Approximations Within the BSE

Once the propagator and kernel are fixed (e.g., at Level 0), there is still a choice of **which excitation channels to include**. This gives two further levels of approximation.

### 3.1 Tamm-Dancoff Approximation (TDA)

The TDA restricts the excitation operator to **forward (particle-hole) processes only**:

```math
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n_0, n} \psi^{n_0 n}_\mathbf{k}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}
```

where $n_0 \in$ occ and $n \in$ unocc. The variational principle leads to a **Hermitian** eigenvalue problem:

```math
\mathcal{A}(\mathbf{q})\, \boldsymbol{\psi}_\mu = \varepsilon_\mu\, \boldsymbol{\psi}_\mu
```

where the $\mathcal{A}$ matrix contains:
- **Diagonal**: mean-field particle-hole energy $E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}$
- **Off-diagonal**: residual interaction (exchange minus direct kernels)

**Properties**:
- Hermitian $\Rightarrow$ real eigenvalues, orthogonal eigenvectors
- All eigenvalues are non-negative (for a stable HF ground state)
- Simple and numerically robust

**Limitation**: ignores ground-state correlations (quantum fluctuations). The HF ground state $|G\rangle$ is treated as the exact vacuum with no zero-point particle-hole pairs.

### 3.2 Random Phase Approximation (RPA)

The RPA extends the excitation operator to include **both forward and backward processes**:

```math
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n_0, n} \left[ X^{n_0 n}_\mathbf{k}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0} - Y^{n_0 n}_\mathbf{k}\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}-\mathbf{q}, n} \right]
```

The backward term $Y$ allows the excitation to **annihilate** virtual particle-hole pairs already present in the correlated ground state. This leads to the **Bosonic BdG** eigenvalue problem:

```math
\begin{pmatrix} \mathcal{A}(\mathbf{q}) & \mathcal{B}(\mathbf{q}) \\ -\mathcal{B}(-\mathbf{q})^* & -\mathcal{A}(-\mathbf{q})^* \end{pmatrix} \begin{pmatrix} X \\ Y \end{pmatrix} = \varepsilon \begin{pmatrix} X \\ Y \end{pmatrix}
```

where:
- $\mathcal{A}(\mathbf{q})$: same TDA matrix (forward-forward coupling)
- $\mathcal{B}(\mathbf{q})$: forward-backward coupling (ground-state correlations)
- The lower-left block $-\mathcal{B}(-\mathbf{q})^*$ and lower-right block $-\mathcal{A}(-\mathbf{q})^*$ follow from the symmetry relations $\mathcal{C}(\mathbf{q}) = -\mathcal{B}(-\mathbf{q})^*$ and $\mathcal{D}(\mathbf{q}) = -\mathcal{A}(-\mathbf{q})^*$, which hold for any Hermitian Hamiltonian

**Properties**:
- Non-Hermitian, but eigenvalues come in $\pm\varepsilon$ pairs
- **Guarantees Goldstone modes**: for spontaneously broken continuous symmetries, the RPA spectrum is exactly gapless at the ordering wavevector. TDA generically gives a small spurious gap.
- Symplectic normalization: $X^\dagger X - Y^\dagger Y = I$

**When RPA matters**:
- Antiferromagnets and other symmetry-broken states (Goldstone theorem)
- Systems where $\mathcal{B}$ is non-negligible (strong ground-state correlations)
- When exact sum rules or conservation laws must be satisfied

### 3.3 Comparison

| Aspect | TDA | RPA |
|:-------|:----|:----|
| Excitation channels | Forward ($ph$) only | Forward ($ph$) + backward ($hp$) |
| Matrix structure | $M \times M$ Hermitian | $2M \times 2M$ non-Hermitian |
| Goldstone theorem | Not guaranteed | Guaranteed |
| Computational cost | $\mathcal{O}(M^3)$ | $\mathcal{O}((2M)^3) \approx 8\times$ TDA |
| Numerical stability | Robust (Hermitian) | Requires care (Cholesky or symplectic diag.) |
| Ground-state correlations | Ignored | Included via $\mathcal{B}$, $\mathcal{D}$ |

### 3.4 Relation to the Single-Mode Approximation

As a historical note, Feynman's **single-mode approximation (SMA)** for superfluid $^4\text{He}$ (1954) can be viewed as a special case of TDA where the envelope function $\psi^n_\mathbf{k}$ is **frozen** to a predetermined form (the density operator $\rho_\mathbf{q}$), leaving no variational freedom. The SMA gives one energy per wavevector via the celebrated formula $\varepsilon(\mathbf{q}) = f(\mathbf{q}) / S(\mathbf{q})$. The TDA/BSE generalizes this by variationally optimizing $\psi$, yielding a complete multi-branch excitation spectrum.

---

## 4. Implementation in MeanFieldTheories.jl

The function `solve_ph_excitations` implements both TDA and RPA at the HF + bare interaction level (Level 0):

- `solver = :TDA`: diagonalizes the Hermitian $\mathcal{A}$ matrix
- `solver = :RPA`: constructs the full Bosonic BdG matrix and solves via Cholesky decomposition (with fallback to direct diagonalization)

The detailed derivations of the $\mathcal{A}$, $\mathcal{B}$, and $\mathcal{D}$ matrices are given in:
- [Tamm-Dancoff Approximation](particle_hole.md): derivation of the $\mathcal{A}$ matrix
- [Random Phase Approximation](particle_hole2.md): derivation of the $\mathcal{B}$ and $\mathcal{D}$ matrices, symplectic normalization, and the Bosonic BdG eigenvalue problem

---

## References

[1] E. E. Salpeter and H. A. Bethe, [A Relativistic Equation for Bound-State Problems](https://doi.org/10.1103/PhysRev.84.1232), Phys. Rev. 84, 1232 (1951).

[2] D. J. Rowe, [Methods for Calculating Ground-State Correlations of Vibrational Nuclei](https://doi.org/10.1103/PhysRev.175.1283), Phys. Rev. 175, 1283 (1968).

[3] G. Onida, L. Reining, and A. Rubio, [Electronic excitations: density-functional versus many-body Green's-function approaches](https://doi.org/10.1103/RevModPhys.74.601), Rev. Mod. Phys. 74, 601 (2002).
