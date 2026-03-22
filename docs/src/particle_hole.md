# Particle-Hole Excitation Theory on Top of the Hartree-Fock Mean Field

---

## 1. Starting Point: The Many-Body Hamiltonian

Consider a general many-body Hamiltonian in second-quantized form:

```math
H = H_0 + H_{\text{int}} = \sum_{ij,ab} t^{ab}_{ij}\,c^\dagger_{ia}\,c_{jb}
    + \sum_{ijkl,abcd} V^{abcd}_{ijkl}\,c^\dagger_{ia}\,c_{jb}\,c^\dagger_{kc}\,c_{ld}
```

where $i,j,k,l$ are site indices and $a,b,c,d$ are orbital (including spin) indices. $t^{ab}_{ij}$ are hopping matrix elements and $V^{abcd}_{ijkl}$ are two-body interaction matrix elements.

After Fourier transformation, the hopping term becomes

```math
H_0 = \sum_{\mathbf{k}, ab} h^{ab}(\mathbf{k})\, c^\dagger_{\mathbf{k}a}\, c_{\mathbf{k}b}
```

where $h^{ab}(\mathbf{k}) = \sum_\delta t^{ab}_\delta\, e^{i\mathbf{k}\cdot\boldsymbol{\delta}}$ ($\delta$ denotes lattice displacement vectors).

The interaction term in momentum space takes the **three-momentum form**:

```math
H_{\text{int}} = \frac{1}{N}\sum_{\mathbf{k}_1\mathbf{k}_2\mathbf{k}_3}\sum_{abcd}
\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)\,
c^\dagger_{\mathbf{k}_1 a}\,c_{\mathbf{k}_2 b}\,c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_4 d}
```

where $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2$ is fixed by momentum conservation, leaving only three independent momenta. $\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)$ is the three-momentum Fourier transform of the interaction potential. The index-operator correspondence is:

- Index $a$: $c^\dagger_{\mathbf{k}_1 a}$
- Index $b$: $c_{\mathbf{k}_2 b}$
- Index $c$: $c^\dagger_{\mathbf{k}_3 c}$
- Index $d$: $c_{\mathbf{k}_4 d}$

---

## 2. Hartree-Fock Mean-Field Approximation

### 2.1 Density Matrix

For a ground state that preserves discrete translational symmetry, the single-particle density matrix is diagonal in momentum space:

```math
\langle c^\dagger_{\mathbf{k}a}\,c_{\mathbf{k}'b}\rangle = \delta_{\mathbf{k},\mathbf{k}'}\,G^{ab}(\mathbf{k})
```

### 2.2 Wick Decomposition

Applying Wick's theorem to the four-fermion operator product (dropping the fully contracted constant), there are four single-contraction channels:

```math
\begin{aligned}
c^\dagger_{\mathbf{k}_1 a}\,c_{\mathbf{k}_2 b}\,c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_4 d}
\;\approx\;&
\underbrace{
+\langle c^\dagger_{\mathbf{k}_1 a} c_{\mathbf{k}_2 b}\rangle\, c^\dagger_{\mathbf{k}_3 c} c_{\mathbf{k}_4 d}
+\langle c^\dagger_{\mathbf{k}_3 c} c_{\mathbf{k}_4 d}\rangle\, c^\dagger_{\mathbf{k}_1 a} c_{\mathbf{k}_2 b}
}_{\text{Hartree (direct)}} \\
&\underbrace{
-\langle c^\dagger_{\mathbf{k}_1 a} c_{\mathbf{k}_4 d}\rangle\, c^\dagger_{\mathbf{k}_3 c} c_{\mathbf{k}_2 b}
-\langle c^\dagger_{\mathbf{k}_3 c} c_{\mathbf{k}_2 b}\rangle\, c^\dagger_{\mathbf{k}_1 a} c_{\mathbf{k}_4 d}
}_{\text{Fock (exchange)}}
\end{aligned}
```

The Hartree terms contract "same-side" operator pairs (the $1,2$ side or the $3,4$ side), while the Fock terms contract "cross-side" pairs, picking up a minus sign from fermionic anticommutation.

### 2.3 Hartree-Fock Self-Energy

Substituting all contractions, the mean-field interaction takes the form

```math
H_{\text{MF}} = \sum_{\mathbf{q},ab} \Sigma^{ab}(\mathbf{q})\, c^\dagger_{\mathbf{q}a}\,c_{\mathbf{q}b}
```

The Hartree-Fock self-energy is

```math
\Sigma^{ab}(\mathbf{q}) = \frac{1}{N}\sum_{\mathbf{k}}\sum_{cd}
\Bigl[
\underbrace{
\widetilde{V}^{cdab}(\mathbf{k},\mathbf{k},\mathbf{q})
+\widetilde{V}^{abcd}(\mathbf{q},\mathbf{q},\mathbf{k})
}_{\text{Hartree}}
\underbrace{
-\widetilde{V}^{cbad}(\mathbf{k},\mathbf{q},\mathbf{q})
-\widetilde{V}^{adcb}(\mathbf{q},\mathbf{k},\mathbf{k})
}_{\text{Fock}}
\Bigr]\,G^{cd}(\mathbf{k})
```

### 2.4 Diagonalization and Band Structure

Define the effective single-particle Hamiltonian matrix

```math
\mathcal{H}^{ab}(\mathbf{k}) = h^{ab}(\mathbf{k}) + \Sigma^{ab}(\mathbf{k})
```

Diagonalize it:

```math
\sum_b \mathcal{H}^{ab}(\mathbf{k})\, U_{bn}(\mathbf{k}) = E^n_\mathbf{k}\, U_{an}(\mathbf{k})
```

where $n$ is the band index, $E^n_\mathbf{k}$ is the mean-field band dispersion, and $U(\mathbf{k})$ is the unitary transformation matrix. Introduce quasiparticle operators:

```math
f_{\mathbf{k}n} = \sum_a U^*_{an}(\mathbf{k})\, c_{\mathbf{k}a}, \qquad c_{\mathbf{k}a} = \sum_n U_{an}(\mathbf{k})\, f_{\mathbf{k}n}
```

The mean-field Hamiltonian in the quasiparticle basis is then diagonal:

```math
\hat{\mathcal{H}}_{\text{MF}} = \sum_{\mathbf{k},n} E^n_\mathbf{k}\, f^\dagger_{\mathbf{k}n}\, f_{\mathbf{k}n} + E_{\text{const}}
```

### 2.5 Hartree-Fock Ground State

The Hartree-Fock ground state $|G\rangle$ is constructed by filling all bands below the Fermi level:

```math
|G\rangle = \prod_{\mathbf{k}} \prod_{n \in \text{occ}} f^\dagger_{\mathbf{k}n}\, |0\rangle
```

where "occ" denotes the set of occupied bands. The density matrix in the quasiparticle basis is

```math
\langle f^\dagger_{\mathbf{k}n}\, f_{\mathbf{k}n'} \rangle = \delta_{nn'}\, \bar{n}_n(\mathbf{k}), \qquad \bar{n}_n(\mathbf{k}) = \begin{cases} 1 & n \in \text{occ} \\ 0 & n \in \text{unocc} \end{cases}
```

---

## 3. Parametrization of the Particle-Hole Excited State

### 3.1 Physical Picture

The simplest excitation above the Hartree-Fock ground state is a **particle-hole pair**: an electron is promoted from an occupied band $n_0$ (the hole) to an unoccupied band $n$ (the particle). In a translationally invariant system, such an excitation carries a definite total momentum $\mathbf{q}$. Depending on the quantum numbers involved, these collective excitations may correspond to magnons, excitons, or other physical objects.

### 3.2 Excited-State Wavefunction

The most general particle-hole excited state at momentum $\mathbf{q}$ is parametrized as

```math
|\mu, \mathbf{q}\rangle = \sum_{\mathbf{k}} \sum_{n_0 \in \text{occ}} \sum_{n \in \text{unocc}} \psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle
```

where:
- $n_0$ runs over all occupied bands (the band where the hole resides)
- $n$ runs over all unoccupied bands (the band where the particle resides)
- $\mathbf{k}$ is the hole momentum (the particle momentum is $\mathbf{k}+\mathbf{q}$)
- $\mu$ is the particle-hole excitation band index (labeling distinct collective modes)
- $\psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})$ is the envelope function, i.e., the variational coefficient to be determined

Define the particle-hole pair operator:

```math
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n_0, n} \psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}
```

so that $|\mu, \mathbf{q}\rangle = \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle$.

### 3.3 Orthonormalization

We require orthonormality among different excited states:

```math
\langle \mu', \mathbf{q} | \mu, \mathbf{q} \rangle = \delta_{\mu\mu'}
```

The inner product is evaluated on the Hartree-Fock ground state. Applying Wick's theorem to the four-operator expectation value

```math
\langle G|\, f^\dagger_{\mathbf{k}', n_0'}\, f_{\mathbf{k}'+\mathbf{q}, n'}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle
```

we use the key facts:
- $n, n'$ are unoccupied bands: $\langle f^\dagger_\cdot\, f_{\mathbf{k}+\mathbf{q}, n} \rangle = 0$
- $n_0, n_0'$ are occupied bands: $\langle f^\dagger_{\mathbf{k}', n_0'}\, f_{\mathbf{k}, n_0} \rangle = \delta_{\mathbf{k}'\mathbf{k}}\,\delta_{n_0' n_0}$

The only surviving fully contracted Wick pairing gives:

```math
\langle G|\, f^\dagger_{\mathbf{k}', n_0'}\, f_{\mathbf{k}'+\mathbf{q}, n'}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle = \delta_{\mathbf{k}'\mathbf{k}}\, \delta_{n_0' n_0}\, \delta_{n'n}
```

Therefore the orthonormalization condition reduces to:

```math
\langle \mu', \mathbf{q} | \mu, \mathbf{q} \rangle = \sum_{\mathbf{k}, n_0, n} \psi^{n_0 n*}_\mathbf{k}(\mu', \mathbf{q})\, \psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q}) = \delta_{\mu\mu'}
```

That is, $\psi^{n_0 n}_\mathbf{k}$ as a vector in the composite $(\mathbf{k}, n_0, n)$ space satisfies the standard orthonormality condition.

---

## 4. Variational Principle and the Eigenvalue Problem

### 4.1 Energy Functional

The excitation energy is computed using the **full original Hamiltonian** $H$ (not the mean-field $H_{\text{MF}}$):

```math
\varepsilon_\mu(\mathbf{q}) = \langle \mu, \mathbf{q} | H | \mu, \mathbf{q} \rangle - E_G
```

where $E_G = \langle G | H | G \rangle$ is the ground-state energy.

### 4.2 Recasting the Excitation Energy in Commutator Form

Using $|\mu, \mathbf{q}\rangle = \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle$, the excitation energy can be written as

```math
\varepsilon_\mu(\mathbf{q}) = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle - E_G
```

To simplify further, we use a simple operator identity. For any operators $A, B$, we have $AB = [A, B] + BA$, and therefore

```math
\hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}} = \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}] + \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H
```

Taking the ground-state expectation value:

```math
\langle G|\, \hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}]\, |G\rangle + \langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H\, |G\rangle
```

Now consider the second term. Within the Hartree-Fock framework, $|G\rangle$ is not the exact eigenstate of $H$, but by the **Brillouin theorem** (Appendix B), $H$ does not mix $|G\rangle$ into any single particle-hole excited state. This means that within the single particle-hole excitation subspace of interest, $H|G\rangle$ effectively behaves as $E_G|G\rangle$. Therefore

```math
\langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H\, |G\rangle = E_G\, \langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = E_G
```

(The last step uses the normalization $\langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = \langle \mu, \mathbf{q}|\mu, \mathbf{q}\rangle = 1$.)

Combining both terms, $E_G$ cancels with $-E_G$, yielding:

```math
\boxed{\varepsilon_\mu(\mathbf{q}) = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}]\, |G\rangle}
```

The advantage of this form is that the commutator $[H, \hat{O}^\dagger]$ is more compact to evaluate than $H\hat{O}^\dagger$ directly, since the commutator naturally reduces the number of operator terms that need to be handled.

### 4.3 Expansion into a Quadratic Form in $\psi^{n_0 n}_\mathbf{k}$

Substituting the definitions of $\hat{O}^\dagger$ and $\hat{O}$ into the above expression. Recall

```math
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{p}, n_0', n'} \psi^{n_0' n'}_\mathbf{p}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}, \qquad \hat{O}_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n_0, n} \psi^{n_0 n*}_\mathbf{k}\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}
```

(The $(\mu, \mathbf{q})$ arguments of $\psi$ are suppressed for brevity.) After substitution, the excitation energy becomes

```math
\varepsilon_\mu(\mathbf{q}) = \sum_{\mathbf{k}, n_0, n}\sum_{\mathbf{p}, n_0', n'} \psi^{n_0 n*}_\mathbf{k}\; \underbrace{\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [H,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]\, |G\rangle}_{\displaystyle \equiv\; \mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})}\; \psi^{n_0' n'}_\mathbf{p}
```

Define the **effective Hamiltonian matrix element**

```math
\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) \equiv \langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [H,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]\, |G\rangle
```

The excitation energy then takes the vector-matrix-vector form:

```math
\varepsilon_\mu(\mathbf{q}) = \sum_{\mathbf{k}, n_0, n}\sum_{\mathbf{p}, n_0', n'} \psi^{n_0 n*}_\mathbf{k}\; \mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\; \psi^{n_0' n'}_\mathbf{p} = \boldsymbol{\psi}^\dagger\, \mathcal{H}^{\text{eff}}(\mathbf{q})\, \boldsymbol{\psi}
```

where $\boldsymbol{\psi}$ is a column vector with composite index $(\mathbf{k}, n_0, n)$. This is a **Hermitian (quadratic) form** in $\psi^{n_0 n}_\mathbf{k}$.

### 4.4 Rayleigh Quotient and the Eigenvalue Problem

Our goal is to find, among all $\boldsymbol{\psi}$ satisfying the normalization constraint $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} = 1$, those that make the excitation energy $\varepsilon$ stationary. This is precisely the standard **Rayleigh quotient** problem.

Define the Rayleigh quotient

```math
R[\boldsymbol{\psi}] = \frac{\boldsymbol{\psi}^\dagger\, \mathcal{H}^{\text{eff}}\, \boldsymbol{\psi}}{\boldsymbol{\psi}^\dagger\, \boldsymbol{\psi}}
```

(Here the denominator involves the $\delta_{\mathbf{k}\mathbf{k}'}\delta_{n_0 n_0'}\delta_{nn'}$ metric, i.e., the standard inner product — this is precisely the significance of the orthonormalization result proved in §3.3. If the metric were not the identity, the problem would become a generalized eigenvalue problem.)

Taking the variation of $R[\boldsymbol{\psi}]$ with respect to $\psi^{n_0 n*}_\mathbf{k}$ and setting $\delta R / \delta \psi^{n_0 n*}_\mathbf{k} = 0$, using the quotient rule:

```math
\frac{\delta}{\delta \psi^{n_0 n*}_\mathbf{k}} \left( \frac{\boldsymbol{\psi}^\dagger \mathcal{H}^{\text{eff}} \boldsymbol{\psi}}{\boldsymbol{\psi}^\dagger \boldsymbol{\psi}} \right) = 0
```

```math
\Longrightarrow \quad \frac{1}{\boldsymbol{\psi}^\dagger \boldsymbol{\psi}} \left[ \sum_{\mathbf{p}, n_0', n'} \mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}\, \psi^{n_0' n'}_\mathbf{p} - R[\boldsymbol{\psi}] \cdot \psi^{n_0 n}_\mathbf{k} \right] = 0
```

Since $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} \neq 0$, the expression in brackets must vanish. Denoting the value of the Rayleigh quotient at the stationary point by $\varepsilon_\mu$, we obtain the standard **eigenvalue problem**:

```math
\boxed{\sum_{\mathbf{p}, n_0', n'} \mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\, \psi^{n_0' n'}_\mathbf{p}(\mu, \mathbf{q}) = \varepsilon_\mu(\mathbf{q})\, \psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})}
```

This is the **Bethe-Salpeter equation** in the Tamm-Dancoff approximation. Its physical content is:

- Diagonalizing the matrix $\mathcal{H}^{\text{eff}}(\mathbf{q})$ yields eigenvalues $\varepsilon_\mu(\mathbf{q})$, the excitation energies of each collective mode $\mu$ at momentum $\mathbf{q}$
- The corresponding eigenvectors $\psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})$ give the particle-hole pair wavefunction of that mode
- The **minimum** of the Rayleigh quotient gives the lowest-energy collective mode; **all stationary points** (saddle points) give the complete excitation spectrum

---

## 5. Complete Derivation of the Effective Hamiltonian $\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})$

### 5.1 Decomposition Strategy

```math
\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = (\mathcal{H}_0)^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + (\mathcal{H}_{\text{int}})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
```

We compute the contributions from the one-body part $[H_0, \cdot]$ and the two-body part $[H_{\text{int}}, \cdot]$ separately.

### 5.2 One-Body Part: $[H_0, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$

Write $H_0$ in the quasiparticle basis:

```math
H_0 = \sum_{\mathbf{p}', mm'} \tilde{h}_{mm'}(\mathbf{p}')\, f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}'m'}
```

where $\tilde{h}_{mm'}(\mathbf{p}') = \sum_{ab} U^*_{am}(\mathbf{p}')\, h^{ab}(\mathbf{p}')\, U_{bm'}(\mathbf{p}')$.

Using the commutation relation for fermionic bilinears (Appendix A):

```math
[f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}'m'},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}] = \delta_{\mathbf{p}', \mathbf{p}+\mathbf{q}}\,\delta_{m'n'}\, f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}, n_0'} - \delta_{\mathbf{p}', \mathbf{p}}\,\delta_{m, n_0'}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}'m'}
```

Summing over $\mathbf{p}', m, m'$:

```math
[H_0,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}] = \sum_m \tilde{h}_{mn'}(\mathbf{p}+\mathbf{q})\, f^\dagger_{\mathbf{p}+\mathbf{q}, m}\, f_{\mathbf{p}, n_0'} - \sum_{m'} \tilde{h}_{n_0' m'}(\mathbf{p})\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, m'}
```

**Physical interpretation**: the first term scatters the particle ($n' \to m$); the second term scatters the hole ($n_0' \to m'$).

Taking the ground-state expectation value $\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n} \cdots |G\rangle$. In the first term, the bra's $f^\dagger_{\mathbf{k},n_0}$ must contract with $f_{\mathbf{p},n_0'}$, requiring $\mathbf{k}=\mathbf{p}$ and $n_0 = n_0'$; then $f_{\mathbf{k}+\mathbf{q},n}$ contracts with $f^\dagger_{\mathbf{p}+\mathbf{q},m}$, giving $m=n$. In the second term, $f_{\mathbf{k}+\mathbf{q},n}$ contracts with $f^\dagger_{\mathbf{p}+\mathbf{q},n'}$, requiring $\mathbf{k}=\mathbf{p}$ and $n=n'$; then $f^\dagger_{\mathbf{k},n_0}$ contracts with $f_{\mathbf{p},m'}$, giving $m' = n_0$. The two terms respectively give:

```math
(\mathcal{H}_0)^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}} \left[ \delta_{n_0 n_0'}\,\tilde{h}_{nn'}(\mathbf{k}+\mathbf{q}) - \delta_{nn'}\, \tilde{h}_{n_0' n_0}(\mathbf{k}) \right]
```

> **Note**: Here $\tilde{h}$ is the bare hopping Hamiltonian in the quasiparticle basis, which is generally **not diagonal**. The full mean-field Hamiltonian $\mathcal{H} = h + \Sigma$ is what diagonalizes. After combining the contributions from $H_0$ and $H_{\text{int}}$ (see §5.5), the mean-field part will fully diagonalize. The hole-scattering term $\tilde{h}_{n_0 n_0'}(\mathbf{k})$ naturally couples different hole bands.

### 5.3 Two-Body Part: Expanding the Commutator

Write $H_{\text{int}}$ in the quasiparticle basis, transforming orbital indices to band indices with the $U$ matrix:

```math
H_{\text{int}} = \frac{1}{N}\sum_{\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}\sum_{m_1 m_2 m_3 m_4} \widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)\, f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}
```

where $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2$, and the interaction in the quasiparticle basis is

```math
\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{p}_1)\, U_{bm_2}(\mathbf{p}_2)\, U^*_{cm_3}(\mathbf{p}_3)\, U_{dm_4}(\mathbf{p}_4)\; \widetilde{V}^{abcd}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)
```

To compute $[H_{\text{int}},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$, use

```math
[f^\dagger_1 f_2 f^\dagger_3 f_4,\; f^\dagger_5 f_6] = [f^\dagger_1 f_2,\; f^\dagger_5 f_6]\, f^\dagger_3 f_4 + f^\dagger_1 f_2\, [f^\dagger_3 f_4,\; f^\dagger_5 f_6]
```

Each bilinear commutator produces two terms (Appendix A), yielding a total of **4 terms**, each containing four fermion operators:

**Part I**: $[f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$

- Term (I-a): $\delta_{\mathbf{p}_2, \mathbf{p}+\mathbf{q}}\,\delta_{m_2,n'}$, producing $f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}, n_0'}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}$
- Term (I-b): $\delta_{\mathbf{p}_1, \mathbf{p}}\,\delta_{m_1,n_0'}$, producing $-f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}$

**Part II**: $[f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$

- Term (II-a): $\delta_{\mathbf{p}_4, \mathbf{p}+\mathbf{q}}\,\delta_{m_4,n'}$, producing $f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}, n_0'}$
- Term (II-b): $\delta_{\mathbf{p}_3, \mathbf{p}}\,\delta_{m_3,n_0'}$, producing $-f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}_4,m_4}$

### 5.4 Ground-State Expectation Values

For each of the four terms, we need to evaluate six-operator expectation values of the form

```math
\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [\text{4 } f/f^\dagger \text{ operators}]\; |G\rangle
```

On the Slater-determinant ground state, Wick's theorem decomposes these into all possible full contractions (three pairs), with the contraction rule

```math
\langle G|\, f^\dagger_{\mathbf{p}m}\, f_{\mathbf{p}'m'}\, |G\rangle = \delta_{\mathbf{p}\mathbf{p}'}\,\delta_{mm'}\, \bar{n}_m(\mathbf{p})
```

**Key selection rules**:

- The external $f_{\mathbf{k}+\mathbf{q}, n}$ ($n \in$ unocc): $\bar{n}_n = 0$, so it cannot pair with any $f^\dagger$ through a Wick contraction — it must "annihilate" against an $f^\dagger$ with matching momentum and band index via the anticommutation relation
- The external $f^\dagger_{\mathbf{k}, n_0}$ ($n_0 \in$ occ): $\bar{n}_{n_0} = 1$, so it can pair normally with a matching $f$

Therefore, in the six-operator expectation value, $f_{\mathbf{k}+\mathbf{q}, n}$ must pair with one of the internal $f^\dagger$ operators, and $f^\dagger_{\mathbf{k}, n_0}$ must pair with one of the internal $f$ operators. The remaining internal pair then contracts with the Fermi sea (self-energy type) or is directly fixed by external quantum numbers (genuine particle-hole scattering).

**Analysis of the four terms** (Term (I-a) detailed; the others below):

#### Term (I-a): Fix $\mathbf{p}_2 = \mathbf{p}+\mathbf{q}$, $m_2 = n'$

Six-operator structure (labeling operator positions 1–6 for sign tracking):

```math
\underset{(1)}{f^\dagger_{\mathbf{k}, n_0}}\, \underset{(2)}{f_{\mathbf{k}+\mathbf{q}, n}}\, \underset{(3)}{f^\dagger_{\mathbf{p}_1,m_1}}\, \underset{(4)}{f_{\mathbf{p}, n_0'}}\, \underset{(5)}{f^\dagger_{\mathbf{p}_3,m_3}}\, \underset{(6)}{f_{\mathbf{p}_4,m_4}}
```

The external annihilator $f_{(2)}$ ($n \in$ unocc) can pair with either of the two internal $f^\dagger$'s: $f^\dagger_{(3)}$ or $f^\dagger_{(5)}$. For each such pairing, the external creator $f^\dagger_{(1)}$ ($n_0 \in$ occ) can pair with any remaining internal $f$. This generates **four** distinct contractions, falling into three classes:

**Contraction (α)**: $(2 \leftrightarrow 3),\; (1 \leftrightarrow 4),\; (5 \leftrightarrow 6)$

Requires $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$; $\mathbf{k} = \mathbf{p}$, $n_0 = n_0'$; $\mathbf{p}_3 = \mathbf{p}_4$, $m_3 = m_4$, $\bar{n}_{m_3} \neq 0$. The sum over $m_3 \in$ occ produces a density matrix — this is part of the Hartree-Fock self-energy correction to the particle energy. [Self-energy type]

**Contraction (β)**: $(2 \leftrightarrow 5),\; (1 \leftrightarrow 4),\; (3 \leftrightarrow 6)$

Requires $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $m_3 = n$; $\mathbf{k} = \mathbf{p}$, $n_0 = n_0'$; $\mathbf{p}_1 = \mathbf{p}_4$, $m_1 = m_4$, $\bar{n}_{m_1} \neq 0$. Also a self-energy correction type. [Self-energy type]

**Contraction (γ-A)**: $(2 \leftrightarrow 3),\; (1 \leftrightarrow 6),\; (5 \leftrightarrow 4)$

Requires $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$; $\mathbf{p}_4 = \mathbf{k}$, $m_4 = n_0$; $\mathbf{p}_3 = \mathbf{p}$, $m_3 = n_0'$, $\bar{n}_{n_0'} = 1$ ✓.

This is a **genuine particle-hole scattering** term — all four band indices ($n, n_0, n', n_0'$) are explicitly fixed, with no Fermi-sea summation.

Momentum check: $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2 = (\mathbf{k}+\mathbf{q}) + \mathbf{p} - (\mathbf{p}+\mathbf{q}) = \mathbf{k}$ ✓.

$V$ tensor: $\widetilde{V}_{m_1, m_2, m_3, m_4} = \widetilde{V}_{n, n', n_0', n_0}(\mathbf{k}+\mathbf{q},\, \mathbf{p}+\mathbf{q},\, \mathbf{p})$.

Sign computation: the pairing permutation $(2,3)(5,4)(1,6)$ relative to the canonical order $(1,2,3,4,5,6)$ has signature $(-1)^3 = -1$.

**Contraction (γ-B)**: $(2 \leftrightarrow 5),\; (1 \leftrightarrow 6),\; (3 \leftrightarrow 4)$

Requires $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $m_3 = n$; $\mathbf{p}_4 = \mathbf{k}$, $m_4 = n_0$; $\mathbf{p}_1 = \mathbf{p}$, $m_1 = n_0'$, $\bar{n}_{n_0'} = 1$ ✓.

Also a **genuine particle-hole scattering** term.

Momentum check: $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2 = \mathbf{p} + (\mathbf{k}+\mathbf{q}) - (\mathbf{p}+\mathbf{q}) = \mathbf{k}$ ✓.

$V$ tensor: $\widetilde{V}_{m_1, m_2, m_3, m_4} = \widetilde{V}_{n_0', n', n, n_0}(\mathbf{p},\, \mathbf{p}+\mathbf{q},\, \mathbf{k}+\mathbf{q})$.

Sign computation: the pairing permutation $(2,5)(3,4)(1,6)$ has signature $(-1)^2 = +1$.

> **Critical observation**: Term (I-a) produces **two** distinct genuine particle-hole scattering contractions (γ-A and γ-B), arising because $f_{(2)}$ can pair with *either* internal $f^\dagger$. These two contractions involve $\widetilde{V}$ evaluated at **different momentum arguments**: $(\mathbf{k}+\mathbf{q},\, \mathbf{p}+\mathbf{q},\, \mathbf{p})$ for γ-A versus $(\mathbf{p},\, \mathbf{p}+\mathbf{q},\, \mathbf{k}+\mathbf{q})$ for γ-B.

#### Term (II-a): Fix $\mathbf{p}_4 = \mathbf{p}+\mathbf{q}$, $m_4 = n'$

Six-operator structure:

```math
\underset{(1)}{f^\dagger_{\mathbf{k}, n_0}}\, \underset{(2)}{f_{\mathbf{k}+\mathbf{q}, n}}\, \underset{(3)}{f^\dagger_{\mathbf{p}_1,m_1}}\, \underset{(4)}{f_{\mathbf{p}_2,m_2}}\, \underset{(5)}{f^\dagger_{\mathbf{p}_3,m_3}}\, \underset{(6)}{f_{\mathbf{p}, n_0'}}
```

By the same logic, there are two self-energy contractions and two γ-type contractions:

**Contraction (γ-A)**: $(2 \leftrightarrow 3),\; (1 \leftrightarrow 4),\; (5 \leftrightarrow 6)$

$\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$; $\mathbf{p}_2 = \mathbf{k}$, $m_2 = n_0$; $\mathbf{p}_3 = \mathbf{p}$, $m_3 = n_0'$, $\bar{n}_{n_0'} = 1$ ✓.

$\mathbf{p}_4 = (\mathbf{k}+\mathbf{q}) + \mathbf{p} - \mathbf{k} = \mathbf{p}+\mathbf{q}$ ✓.

$V$ tensor: $\widetilde{V}_{n, n_0, n_0', n'}(\mathbf{k}+\mathbf{q},\, \mathbf{k},\, \mathbf{p})$.

Sign: the pairing permutation $(2,3)(1,4)(5,6)$ has signature $(-1)^2 = +1$.

**Contraction (γ-B)**: $(2 \leftrightarrow 5),\; (1 \leftrightarrow 4),\; (3 \leftrightarrow 6)$

$\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $m_3 = n$; $\mathbf{p}_2 = \mathbf{k}$, $m_2 = n_0$; $\mathbf{p}_1 = \mathbf{p}$, $m_1 = n_0'$, $\bar{n}_{n_0'} = 1$ ✓.

$\mathbf{p}_4 = \mathbf{p} + (\mathbf{k}+\mathbf{q}) - \mathbf{k} = \mathbf{p}+\mathbf{q}$ ✓.

$V$ tensor: $\widetilde{V}_{n_0', n_0, n, n'}(\mathbf{p},\, \mathbf{k},\, \mathbf{k}+\mathbf{q})$.

Sign: the pairing permutation $(2,5)(1,4)(3,6)$ has signature $(-1)^3 = -1$.

#### Terms (I-b) and (II-b): No genuine particle-hole scattering contributions

In Term (I-b), the delta fixes $\mathbf{p}_1 = \mathbf{p}$, $m_1 = n_0'$, producing the internal operator string $-f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}$. Any γ-type contraction that does not fix $\mathbf{k} = \mathbf{p}$ necessarily leaves the pair $f^\dagger_{\mathbf{p}+\mathbf{q}, n'}$ to contract with the Fermi sea. Since $n' \in$ unocc implies $\bar{n}_{n'} = 0$, these contractions vanish identically.

The same argument applies to Term (II-b), where the delta fixes $\mathbf{p}_3 = \mathbf{p}$, $m_3 = n_0'$, leaving $f^\dagger_{\mathbf{p}+\mathbf{q}, n'}$ (unocc) in the remaining internal pair.

**Therefore, the genuine particle-hole scattering kernel receives contributions only from Terms (I-a) and (II-a), each of which contributes two γ contractions.**

### 5.5 Combining the Results

**(A) Mean-field part**: From contractions of type (α) and (β) (one external operator pair matched + one internal pair contracted with the Fermi sea). These require $\mathbf{k} = \mathbf{p}$ and $n_0 = n_0'$. They exactly reproduce the Hartree-Fock self-energy corrections to the particle and hole. Combined with the one-body part from §5.2, $\tilde{h}$ is replaced by the full mean-field Hamiltonian $\tilde{\mathcal{H}} = \tilde{h} + \tilde{\Sigma}$, which fully diagonalizes to the mean-field energies:

```math
(\mathcal{H}_{\text{MF}})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right)
```

This is the free particle-hole pair energy: particle energy minus hole energy.

**(B) Particle-hole interaction kernel**: From the four γ contractions identified in §5.4 — two from Term (I-a) and two from Term (II-a). These involve no Fermi-sea summation and all four band indices ($n, n_0, n', n_0'$) are explicitly fixed. **These terms do not require $n_0 = n_0'$** — the interaction kernel couples different hole bands.

### 5.6 Two Topologies of the Particle-Hole Interaction Kernel

The four γ contractions are classified into two topologies based on the **momentum transfer within each bilinear pair** of the interaction vertex $c^\dagger_{\mathbf{k}_1 a}\, c_{\mathbf{k}_2 b}\, c^\dagger_{\mathbf{k}_3 c}\, c_{\mathbf{k}_4 d}$:

- **Direct (density) topology**: momentum transfer $\mathbf{k}_1 - \mathbf{k}_2 = \mathbf{k} - \mathbf{p}$, the relative momentum between hole and particle. The particle and hole lines run in parallel.
- **Exchange topology**: momentum transfer $\mathbf{k}_1 - \mathbf{k}_2 = \pm\mathbf{q}$, the total particle-hole pair momentum. The particle and hole lines cross.

Inspecting the $V$ arguments of each γ contraction:

| Contraction | $V$ momenta $(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)$ | $\mathbf{k}_1 - \mathbf{k}_2$ | Topology | Sign |
|---|---|---|---|---|
| (I-a)-γA | $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$ | $\mathbf{k}-\mathbf{p}$ | Direct | $-1$ |
| (II-a)-γB | $(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$ | $\mathbf{p}-\mathbf{k}$ | Direct | $-1$ |
| (II-a)-γA | $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$ | $\mathbf{q}$ | Exchange | $+1$ |
| (I-a)-γB | $(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$ | $-\mathbf{q}$ | Exchange | $+1$ |

Each channel receives **two contributions** involving $\widetilde{V}$ at **different momentum arguments**.

**(B1) Direct channel (density channel)** $\mathcal{K}^{\text{d}}$:

From (I-a)-γA and (II-a)-γB. Physical process: the interaction line connects the particle and hole lines, each scattering without exchange.

```math
\boxed{
\begin{aligned}
\mathcal{K}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = -\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{bn'}(\mathbf{p}+\mathbf{q})\, U^*_{cn_0'}(\mathbf{p})\, U_{dn_0}(\mathbf{k})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q},\, \mathbf{p}+\mathbf{q},\, \mathbf{p}) \\
+\; & U^*_{an_0'}(\mathbf{p})\, U_{bn_0}(\mathbf{k})\, U^*_{cn}(\mathbf{k}+\mathbf{q})\, U_{dn'}(\mathbf{p}+\mathbf{q})\; \widetilde{V}^{abcd}(\mathbf{p},\, \mathbf{k},\, \mathbf{k}+\mathbf{q})
\Big]
\end{aligned}
}
```

Feynman diagram (both contributions have the same topology):

```
  particle: (p+q, n') ──[V]──▸ (k+q, n)
                            |
  hole:     (k, n₀)   ──[V]──▸ (p, n₀')
```

**First term** (from (I-a)-γA): $\widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q}, \mathbf{p}+\mathbf{q}, \mathbf{p})$, where the first bilinear pair $(a,b)$ of the interaction connects particle-out and particle-in, and the second pair $(c,d)$ connects hole-out and hole-in.

**Second term** (from (II-a)-γB): $\widetilde{V}^{abcd}(\mathbf{p}, \mathbf{k}, \mathbf{k}+\mathbf{q})$, where the first pair $(a,b)$ connects hole-out and hole-in, and the second pair $(c,d)$ connects particle-out and particle-in. This is the same physical process with the roles of the two bilinear pairs in $H_{\text{int}}$ swapped.

Momentum conservation checks:

- First term: $\mathbf{k}_1 = \mathbf{k}+\mathbf{q}$, $\mathbf{k}_2 = \mathbf{p}+\mathbf{q}$, $\mathbf{k}_3 = \mathbf{p}$, $\mathbf{k}_4 = \mathbf{k}$. ✓
- Second term: $\mathbf{k}_1 = \mathbf{p}$, $\mathbf{k}_2 = \mathbf{k}$, $\mathbf{k}_3 = \mathbf{k}+\mathbf{q}$, $\mathbf{k}_4 = \mathbf{p}+\mathbf{q}$. ✓

**(B2) Exchange channel** $\mathcal{K}^{\text{x}}$:

From (II-a)-γA and (I-a)-γB. Physical process: the particle and hole lines cross through the interaction vertex.

```math
\boxed{
\begin{aligned}
\mathcal{K}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = +\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{bn_0}(\mathbf{k})\, U^*_{cn_0'}(\mathbf{p})\, U_{dn'}(\mathbf{p}+\mathbf{q})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q},\, \mathbf{k},\, \mathbf{p}) \\
+\; & U^*_{an_0'}(\mathbf{p})\, U_{bn'}(\mathbf{p}+\mathbf{q})\, U^*_{cn}(\mathbf{k}+\mathbf{q})\, U_{dn_0}(\mathbf{k})\; \widetilde{V}^{abcd}(\mathbf{p},\, \mathbf{p}+\mathbf{q},\, \mathbf{k}+\mathbf{q})
\Big]
\end{aligned}
}
```

Feynman diagram (both contributions have the same topology):

```
  particle: (p+q, n') ──╲    ╱──▸ (k+q, n)
                          ╲  ╱
                           ╳  [V]
                          ╱  ╲
  hole:     (k, n₀)   ──╱    ╲──▸ (p, n₀')
```

**First term** (from (II-a)-γA): $\widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q}, \mathbf{k}, \mathbf{p})$, where the first pair $(a,b)$ connects particle-out to hole-in, and the second pair $(c,d)$ connects hole-out to particle-in.

**Second term** (from (I-a)-γB): $\widetilde{V}^{abcd}(\mathbf{p}, \mathbf{p}+\mathbf{q}, \mathbf{k}+\mathbf{q})$, the same topology with swapped bilinear pair roles.

Momentum conservation checks:

- First term: $\mathbf{k}_1 = \mathbf{k}+\mathbf{q}$, $\mathbf{k}_2 = \mathbf{k}$, $\mathbf{k}_3 = \mathbf{p}$, $\mathbf{k}_4 = \mathbf{p}+\mathbf{q}$. ✓
- Second term: $\mathbf{k}_1 = \mathbf{p}$, $\mathbf{k}_2 = \mathbf{p}+\mathbf{q}$, $\mathbf{k}_3 = \mathbf{k}+\mathbf{q}$, $\mathbf{k}_4 = \mathbf{k}$. ✓

> **Note on signs**: The direct channel carries a minus sign and the exchange channel carries a plus sign, arising from the relative fermion signs in the different contraction topologies. The specific signs depend on the operator ordering convention ($c^\dagger c c^\dagger c$) of the original Hamiltonian and the number of fermion line crossings in the Wick contractions.

> **Physical origin of the two terms per channel**: The Hamiltonian $H_{\text{int}} = \frac{1}{N}\sum V^{abcd} c^\dagger_a c_b c^\dagger_c c_d$ contains two bilinear pairs: $(c^\dagger_a c_b)$ and $(c^\dagger_c c_d)$. When computing $[H_{\text{int}}, f^\dagger f]$, Part I of the commutator (§5.3) acts on the first pair, while Part II acts on the second pair. For each topology (direct or exchange), the external operators can match against *either* bilinear pair of $H_{\text{int}}$, generating two contributions with different $V$ momentum arguments. This is a direct consequence of the $c^\dagger c c^\dagger c$ ordering — unlike the normal-ordered $\frac{1}{2}c^\dagger c^\dagger cc$ form where antisymmetry absorbs this doubling into a factor of 2.

### 5.6.1 Simplification for On-Site Interactions

For **on-site interactions** (e.g., Hubbard, Kanamori), $\widetilde{V}^{abcd}$ is independent of momenta. In this case, the second term in each channel can be related to the first by relabeling dummy orbital indices $a \leftrightarrow c$, $b \leftrightarrow d$:

```math
\text{Second term of } \mathcal{K}^{\text{d}}: \quad \sum_{abcd} [\cdots]\, \widetilde{V}^{abcd} \xrightarrow{a\leftrightarrow c,\, b\leftrightarrow d} \sum_{abcd} [\cdots]\, \widetilde{V}^{cdab}
```

so the direct kernel becomes

```math
\mathcal{K}^{\text{d}} \propto \sum_{abcd} U^*_{an}\, U_{bn'}\, U^*_{cn_0'}\, U_{dn_0} \left[\widetilde{V}^{abcd} + \widetilde{V}^{cdab}\right]
```

and similarly for $\mathcal{K}^{\text{x}}$. This is equivalent to replacing $\widetilde{V}$ by a "particle-exchange-symmetrized" interaction $\bar{V}^{abcd} = \widetilde{V}^{abcd} + \widetilde{V}^{cdab}$.

> **Warning**: For **momentum-dependent interactions** (including on-site interactions in a reduced Brillouin zone, e.g., after magnetic unit cell folding), the two $V$ terms have genuinely different momentum arguments and **cannot** be combined by simple relabeling. They must be evaluated separately.

### 5.7 Final Result

```math
\boxed{
\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right) + \mathcal{K}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + \mathcal{K}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
}
```

where:

- **First term**: free particle-hole pair energy (entirely from the mean field), diagonal in all indices
- **Second term**: direct channel kernel $\mathcal{K}^{\text{d}}$, with two $\widetilde{V}$ contributions at momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$ and $(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$
- **Third term**: exchange channel kernel $\mathcal{K}^{\text{x}}$, with two $\widetilde{V}$ contributions at momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$ and $(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$

In the general case both channels contribute. Depending on the symmetries of a specific model, one of them may vanish. **Both kernel terms couple different hole bands** ($n_0 \neq n_0'$), which is essential for correctly describing collective excitations involving multiple valence bands.

---

## 6. Relation to the Tamm-Dancoff and Random Phase Approximations

### 6.1 Tamm-Dancoff Approximation (TDA)

The eigenvalue problem derived above,

```math
\sum_{\mathbf{p}, n_0', n'} \mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\, \psi^{n_0' n'}_\mathbf{p}(\mu, \mathbf{q}) = \varepsilon_\mu(\mathbf{q})\, \psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})
```

which includes only particle-hole ($ph$) excitations, is equivalent to the **Tamm-Dancoff approximation (TDA)**.

### 6.2 Full RPA

The full **random phase approximation (RPA)** additionally includes hole-particle ($hp$) backward excitations:

```math
|\mu, \mathbf{q}\rangle_{\text{RPA}} = \sum_{\mathbf{k}, n_0, n} \left[ X^{n_0 n}_\mathbf{k}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0} - Y^{n_0 n}_\mathbf{k}\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n} \right] |G\rangle
```

This leads to the standard matrix form of the RPA equation:

```math
\begin{pmatrix} A & B \\ -B^* & -A^* \end{pmatrix} \begin{pmatrix} X \\ Y \end{pmatrix} = \varepsilon \begin{pmatrix} X \\ Y \end{pmatrix}
```

where the $A$ matrix is precisely $\mathcal{H}^{\text{eff}}$ (i.e., the TDA matrix), and the $B$ matrix describes ground-state correlations. The Tamm-Dancoff approximation corresponds to $B = 0$.

### 6.3 Collective Excitation Dispersion

The eigenvalues $\varepsilon_\mu(\mathbf{q})$ give the dispersion relation of collective excitations, with $\mu$ labeling different excitation bands. If the system spontaneously breaks a continuous symmetry, the Goldstone theorem guarantees the existence of a gapless collective mode at the ordering wavevector $\mathbf{q} = \mathbf{Q}$.

---

## 7. Summary: Complete Computational Workflow from the Original Hamiltonian to the Excitation Spectrum

1. **Hartree-Fock self-consistent calculation**:
   - Initialize the density matrix $G^{ab}(\mathbf{k})$
   - Construct $\mathcal{H}^{ab}(\mathbf{k}) = h^{ab}(\mathbf{k}) + \Sigma^{ab}(\mathbf{k})$
   - Diagonalize to obtain $E^n_\mathbf{k}$, $U(\mathbf{k})$
   - Update $G^{ab}(\mathbf{k})$, iterate until convergence

2. **Construct the effective Hamiltonian**:
   - For each $\mathbf{q}$, build the matrix $\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})$ with composite index $(\mathbf{k}, n_0, n)$
   - Diagonal part: $\delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k})$
   - Off-diagonal part: $\mathcal{K}^{\text{d}} + \mathcal{K}^{\text{x}}$, each containing two $\widetilde{V}$ terms at different momentum arguments (see §5.6)

3. **Diagonalize the effective Hamiltonian**:
   - $\mathcal{H}^{\text{eff}}\, \psi = \varepsilon\, \psi$
   - Eigenvalues $\varepsilon_\mu(\mathbf{q})$ give the excitation spectrum
   - Eigenvectors $\psi^{n_0 n}_\mathbf{k}(\mu, \mathbf{q})$ give the excited-state wavefunctions

4. **Physical observables**:
   - Dynamic structure factor: $S(\mathbf{q}, \omega) = \sum_\mu |\langle \mu, \mathbf{q} | \hat{A}_\mathbf{q} | G \rangle|^2\, \delta(\omega - \varepsilon_\mu(\mathbf{q}))$
   - where $\hat{A}_\mathbf{q}$ is the operator that couples to the excitation of interest (e.g., spin operator, density operator); its matrix elements are computed via $U$ and $\psi$

---

## Appendix A: Useful Commutation Relations

For fermionic bilinear operators:

```math
[f^\dagger_\alpha f_\beta,\; f^\dagger_\gamma f_\delta] = \delta_{\beta\gamma}\, f^\dagger_\alpha f_\delta - \delta_{\alpha\delta}\, f^\dagger_\gamma f_\beta
```

Derivation: using $\{f_\alpha, f^\dagger_\beta\} = \delta_{\alpha\beta}$,

```math
\begin{aligned}
f^\dagger_\alpha f_\beta f^\dagger_\gamma f_\delta &= f^\dagger_\alpha (\delta_{\beta\gamma} - f^\dagger_\gamma f_\beta) f_\delta = \delta_{\beta\gamma}\, f^\dagger_\alpha f_\delta - f^\dagger_\alpha f^\dagger_\gamma f_\beta f_\delta \\
f^\dagger_\gamma f_\delta f^\dagger_\alpha f_\beta &= \delta_{\delta\alpha}\, f^\dagger_\gamma f_\beta - f^\dagger_\gamma f^\dagger_\alpha f_\delta f_\beta = \delta_{\delta\alpha}\, f^\dagger_\gamma f_\beta - f^\dagger_\alpha f^\dagger_\gamma f_\beta f_\delta
\end{aligned}
```

(The last step uses $f^\dagger_\gamma f^\dagger_\alpha = -f^\dagger_\alpha f^\dagger_\gamma$ and $f_\delta f_\beta = -f_\beta f_\delta$.) Subtracting the two expressions yields the result.

## Appendix B: Brillouin's Theorem

The matrix element of the full Hamiltonian $H$ between the Hartree-Fock ground state $|G\rangle$ and any single particle-hole excited state vanishes:

```math
\langle G | H | f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0} | G \rangle = 0 \quad (n \in \text{unocc}, \; n_0 \in \text{occ})
```

This is a direct consequence of the Hartree-Fock equations (the self-consistency condition). Physically, the Hartree-Fock ground state is the "optimal" single Slater determinant, whose energy cannot be lowered by any single particle-hole excitation.

This theorem guarantees the validity of recasting $\varepsilon$ in commutator form in §4.2 ($H|G\rangle = E_G|G\rangle$ is effective within the single particle-hole subspace).

## Appendix C: $U$-Matrix Transformation and the Orbital-to-Band Conversion of $\widetilde{V}$

The interaction $\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)$ in the original orbital space is transformed to the band basis via the $U$ matrix:

```math
\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{k}_1)\, U_{bm_2}(\mathbf{k}_2)\, U^*_{cm_3}(\mathbf{k}_3)\, U_{dm_4}(\mathbf{k}_4)\; \widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)
```

where $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2$.

The convention is that each orbital index in $\widetilde{V}^{abcd}$ is contracted with the $U$ matrix evaluated at the **corresponding momentum**: $a$ with $U^*(\mathbf{k}_1)$, $b$ with $U(\mathbf{k}_2)$, $c$ with $U^*(\mathbf{k}_3)$, $d$ with $U(\mathbf{k}_4)$. The boxed formulas for $\mathcal{K}^{\text{d}}$ and $\mathcal{K}^{\text{x}}$ in §5.6 follow this convention exactly: each of the two $\widetilde{V}$ terms within each kernel has its own momentum assignment and its own $U$-matrix mapping.

When constructing $\mathcal{K}^{\text{d}}$ and $\mathcal{K}^{\text{x}}$, only the components of $\widetilde{V}$ with specific band indices ($n, n', n_0, n_0'$) are needed. In practice, there is no need to fully transform the entire tensor — one can directly contract with the relevant columns of the $U$ matrix in orbital space. Note that each kernel requires **two calls** to the interaction function $\widetilde{V}^{abcd}(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)$ with different momentum arguments.

**Implementation note**: When implementing the kernels as matrix multiplications, care must be taken with the index grouping. For a given $\widetilde{V}^{abcd}$ term, one must identify which orbital indices are "free" (to be projected onto bands by $U^\dagger$ and $U$) and which are "contracted" (summed against specific columns of $U$). The free-index pair differs between the direct and exchange topologies:

- **Direct channel**: the first pair $(a,b)$ maps to (particle-out, particle-in) or (hole-out, hole-in) — in either case, the free indices for the output matrix are $(a,b)$, and the contracted indices are $(c,d)$.
- **Exchange channel**: the first pair $(a,b)$ maps to (particle-out, hole-in) or (hole-out, particle-in) — the free indices for the output matrix are $(a,d)$, and the contracted indices are $(b,c)$. To use a matrix multiplication, the tensor must be permuted to $(a,d,b,c)$ order before reshaping.

Confusing the free-index grouping between channels is equivalent to applying the direct-channel (bubble) topology to the exchange (crossed) diagram, or vice versa.