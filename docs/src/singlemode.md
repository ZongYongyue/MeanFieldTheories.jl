# Single-Mode Approximation

This page introduces a theory of collective excitations on top of a
Hartree-Fock ground state, which is commonly called the **single-mode
approximation** (SMA) or **Bethe-Salpeter
equation** (BSE).


---

## 1. Starting Point: The Many-Body Hamiltonian

Consider a general many-body Hamiltonian in second-quantized form:

$$
H = H_0 + H_{\text{int}} = \sum_{ij,ab} t^{ab}_{ij}\,c^\dagger_{ia}\,c_{jb}
    + \sum_{ijkl,abcd} V^{abcd}_{ijkl}\,c^\dagger_{ia}\,c_{jb}\,c^\dagger_{kc}\,c_{ld}
$$

where $i,j,k,l$ are site indices and $a,b,c,d$ are orbital (including spin) indices. $t^{ab}_{ij}$ are hopping matrix elements and $V^{abcd}_{ijkl}$ are two-body interaction matrix elements.

After Fourier transformation, the hopping term becomes

$$
H_0 = \sum_{\mathbf{k}, ab} h^{ab}(\mathbf{k})\, c^\dagger_{\mathbf{k}a}\, c_{\mathbf{k}b}
$$

where $h^{ab}(\mathbf{k}) = \sum_\delta t^{ab}_\delta\, e^{i\mathbf{k}\cdot\boldsymbol{\delta}}$ ($\delta$ denotes lattice displacement vectors).

The interaction term in momentum space takes the **three-momentum form**:

$$
H_{\text{int}} = \frac{1}{N}\sum_{\mathbf{k}_1\mathbf{k}_2\mathbf{k}_3}\sum_{abcd}
\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)\,
c^\dagger_{\mathbf{k}_1 a}\,c_{\mathbf{k}_2 b}\,c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_4 d}
$$

where $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2$ is fixed by momentum conservation, leaving only three independent momenta. $\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)$ is the three-momentum Fourier transform of the interaction potential. The index-operator correspondence is:

- Index $a$: $c^\dagger_{\mathbf{k}_1 a}$
- Index $b$: $c_{\mathbf{k}_2 b}$
- Index $c$: $c^\dagger_{\mathbf{k}_3 c}$
- Index $d$: $c_{\mathbf{k}_4 d}$

---

## 2. Hartree-Fock Mean-Field Approximation

### 2.1 Density Matrix

For a ground state that preserves discrete translational symmetry, the single-particle density matrix is diagonal in momentum space:

$$
\langle c^\dagger_{\mathbf{k}a}\,c_{\mathbf{k}'b}\rangle = \delta_{\mathbf{k},\mathbf{k}'}\,G^{ab}(\mathbf{k})
$$

### 2.2 Wick Decomposition

Applying Wick's theorem to the four-fermion operator product (dropping the fully contracted constant), there are four single-contraction channels:

$$
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
$$

The Hartree terms contract "same-side" operator pairs (the $1,2$ side or the $3,4$ side), while the Fock terms contract "cross-side" pairs, picking up a minus sign from fermionic anticommutation.

### 2.3 Hartree-Fock Self-Energy

Substituting all contractions, the mean-field interaction takes the form

$$
H_{\text{MF}} = \sum_{\mathbf{q},ab} \Sigma^{ab}(\mathbf{q})\, c^\dagger_{\mathbf{q}a}\,c_{\mathbf{q}b}
$$

The Hartree-Fock self-energy is

$$
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
$$

### 2.4 Diagonalization and Band Structure

Define the effective single-particle Hamiltonian matrix

$$
\mathcal{H}^{ab}(\mathbf{k}) = h^{ab}(\mathbf{k}) + \Sigma^{ab}(\mathbf{k})
$$

Diagonalize it:

$$
\sum_b \mathcal{H}^{ab}(\mathbf{k})\, U_{bn}(\mathbf{k}) = E^n_\mathbf{k}\, U_{an}(\mathbf{k})
$$

where $n$ is the band index, $E^n_\mathbf{k}$ is the mean-field band dispersion, and $U(\mathbf{k})$ is the unitary transformation matrix. Introduce quasiparticle operators:

$$
f_{\mathbf{k}n} = \sum_a U^*_{an}(\mathbf{k})\, c_{\mathbf{k}a}, \qquad c_{\mathbf{k}a} = \sum_n U_{an}(\mathbf{k})\, f_{\mathbf{k}n}
$$

The mean-field Hamiltonian in the quasiparticle basis is then diagonal:

$$
\hat{\mathcal{H}}_{\text{MF}} = \sum_{\mathbf{k},n} E^n_\mathbf{k}\, f^\dagger_{\mathbf{k}n}\, f_{\mathbf{k}n} + E_{\text{const}}
$$

### 2.5 Hartree-Fock Ground State

The Hartree-Fock ground state $|G\rangle$ is constructed by filling all bands below the Fermi level:

$$
|G\rangle = \prod_{\mathbf{k}} \prod_{n \in \text{occ}} f^\dagger_{\mathbf{k}n}\, |0\rangle
$$

where "occ" denotes the set of occupied bands. The density matrix in the quasiparticle basis is

$$
\langle f^\dagger_{\mathbf{k}n}\, f_{\mathbf{k}n'} \rangle = \delta_{nn'}\, \bar{n}_n(\mathbf{k}), \qquad \bar{n}_n(\mathbf{k}) = \begin{cases} 1 & n \in \text{occ} \\ 0 & n \in \text{unocc} \end{cases}
$$

---

## 3. Parametrization of the Particle-Hole Excited State

### 3.1 Physical Picture

The simplest excitation above the Hartree-Fock ground state is a **particle-hole pair**: an electron is promoted from an occupied band $n_0$ (the hole) to an unoccupied band $n$ (the particle). In a translationally invariant system, such an excitation carries a definite total momentum $\mathbf{q}$. Depending on the quantum numbers involved, these collective excitations may correspond to magnons, excitons, or other physical objects.

### 3.2 Excited-State Wavefunction

The particle-hole excited state is parametrized as

$$
|\mu, \mathbf{q}\rangle = \sum_{\mathbf{k}, n \in \text{unocc}} \psi^n_\mathbf{k}(\mu, \mathbf{q})\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle
$$

where:
- $n_0$ is a fixed occupied band index (the band where the hole resides)
- $n$ runs over all unoccupied bands (the band where the particle resides)
- $\mathbf{k}$ is the hole momentum (the particle momentum is $\mathbf{k}+\mathbf{q}$)
- $\mu$ is the particle-hole excitation band index (labeling distinct collective modes)
- $\psi^n_\mathbf{k}(\mu, \mathbf{q})$ is the envelope function, i.e., the variational coefficient to be determined

Define the particle-hole pair operator:

$$
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n} \psi^n_\mathbf{k}(\mu, \mathbf{q})\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}
$$

so that $|\mu, \mathbf{q}\rangle = \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle$.

### 3.3 Orthonormalization

We require orthonormality among different excited states:

$$
\langle \mu', \mathbf{q} | \mu, \mathbf{q} \rangle = \delta_{\mu\mu'}
$$

The inner product is evaluated on the Hartree-Fock ground state. Applying Wick's theorem to the four-operator expectation value

$$
\langle G|\, f^\dagger_{\mathbf{k}', n_0}\, f_{\mathbf{k}'+\mathbf{q}, n'}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle
$$

we use the key facts:
- $n, n'$ are unoccupied bands: $\langle f^\dagger_\cdot\, f_{\mathbf{k}+\mathbf{q}, n} \rangle = 0$
- $n_0$ is an occupied band: $\langle f^\dagger_{\mathbf{k}', n_0}\, f_{\mathbf{k}, n_0} \rangle = \delta_{\mathbf{k}'\mathbf{k}}$

The only surviving fully contracted Wick pairing gives:

$$
\langle G|\, f^\dagger_{\mathbf{k}', n_0}\, f_{\mathbf{k}'+\mathbf{q}, n'}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0}\, |G\rangle = \delta_{\mathbf{k}'\mathbf{k}}\, \delta_{n'n}
$$

Therefore the orthonormalization condition reduces to:

$$
\langle \mu', \mathbf{q} | \mu, \mathbf{q} \rangle = \sum_{\mathbf{k}, n} \psi^{n*}_\mathbf{k}(\mu', \mathbf{q})\, \psi^n_\mathbf{k}(\mu, \mathbf{q}) = \delta_{\mu\mu'}
$$

That is, $\psi^n_\mathbf{k}$ as a vector in the composite $(\mathbf{k}, n)$ space satisfies the standard orthonormality condition.

---

## 4. Variational Principle and the Eigenvalue Problem

### 4.1 Energy Functional

The excitation energy is computed using the **full original Hamiltonian** $H$ (not the mean-field $H_{\text{MF}}$):

$$
\varepsilon_\mu(\mathbf{q}) = \langle \mu, \mathbf{q} | H | \mu, \mathbf{q} \rangle - E_G
$$

where $E_G = \langle G | H | G \rangle$ is the ground-state energy.

### 4.2 Recasting the Excitation Energy in Commutator Form

Using $|\mu, \mathbf{q}\rangle = \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle$, the excitation energy can be written as

$$
\varepsilon_\mu(\mathbf{q}) = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle - E_G
$$

To simplify further, we use a simple operator identity. For any operators $A, B$, we have $AB = [A, B] + BA$, and therefore

$$
\hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}} = \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}] + \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H
$$

Taking the ground-state expectation value:

$$
\langle G|\, \hat{O}_{\mu\mathbf{q}}\, H\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}]\, |G\rangle + \langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H\, |G\rangle
$$

Now consider the second term. Within the Hartree-Fock framework, $|G\rangle$ is not the exact eigenstate of $H$, but by the **Brillouin theorem** (Appendix B), $H$ does not mix $|G\rangle$ into any single particle-hole excited state. This means that within the single particle-hole excitation subspace of interest, $H|G\rangle$ effectively behaves as $E_G|G\rangle$. Therefore

$$
\langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, H\, |G\rangle = E_G\, \langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = E_G
$$

(The last step uses the normalization $\langle G|\, \hat{O}_{\mu\mathbf{q}}\, \hat{O}^\dagger_{\mu\mathbf{q}}\, |G\rangle = \langle \mu, \mathbf{q}|\mu, \mathbf{q}\rangle = 1$.)

Combining both terms, $E_G$ cancels with $-E_G$, yielding:

$$
\boxed{\varepsilon_\mu(\mathbf{q}) = \langle G|\, \hat{O}_{\mu\mathbf{q}}\, [H, \hat{O}^\dagger_{\mu\mathbf{q}}]\, |G\rangle}
$$

The advantage of this form is that the commutator $[H, \hat{O}^\dagger]$ is more compact to evaluate than $H\hat{O}^\dagger$ directly, since the commutator naturally reduces the number of operator terms that need to be handled.

### 4.3 Expansion into a Quadratic Form in $\psi^n_\mathbf{k}$

Substituting the definitions of $\hat{O}^\dagger$ and $\hat{O}$ into the above expression. Recall

$$
\hat{O}^\dagger_{\mu\mathbf{q}} = \sum_{\mathbf{p}, n'} \psi^{n'}_\mathbf{p}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}, \qquad \hat{O}_{\mu\mathbf{q}} = \sum_{\mathbf{k}, n} \psi^{n*}_\mathbf{k}\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}
$$

(The $(\mu, \mathbf{q})$ arguments of $\psi$ are suppressed for brevity.) After substitution, the excitation energy becomes

$$
\varepsilon_\mu(\mathbf{q}) = \sum_{\mathbf{k}, n}\sum_{\mathbf{p}, n'} \psi^{n*}_\mathbf{k}\; \underbrace{\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [H,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]\, |G\rangle}_{\displaystyle \equiv\; \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})}\; \psi^{n'}_\mathbf{p}
$$

Define the **effective Hamiltonian matrix element**

$$
\mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) \equiv \langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [H,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]\, |G\rangle
$$

The excitation energy then takes the vector-matrix-vector form:

$$
\varepsilon_\mu(\mathbf{q}) = \sum_{\mathbf{k}, n}\sum_{\mathbf{p}, n'} \psi^{n*}_\mathbf{k}\; \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\; \psi^{n'}_\mathbf{p} = \boldsymbol{\psi}^\dagger\, \mathcal{H}^{\text{eff}}(\mathbf{q})\, \boldsymbol{\psi}
$$

where $\boldsymbol{\psi}$ is a column vector with composite index $(\mathbf{k}, n)$. This is a **Hermitian (quadratic) form** in $\psi^n_\mathbf{k}$.

### 4.4 Rayleigh Quotient and the Eigenvalue Problem

Our goal is to find, among all $\boldsymbol{\psi}$ satisfying the normalization constraint $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} = 1$, those that make the excitation energy $\varepsilon$ stationary. This is precisely the standard **Rayleigh quotient** problem.

Define the Rayleigh quotient

$$
R[\boldsymbol{\psi}] = \frac{\boldsymbol{\psi}^\dagger\, \mathcal{H}^{\text{eff}}\, \boldsymbol{\psi}}{\boldsymbol{\psi}^\dagger\, \boldsymbol{\psi}}
$$

(Here the denominator involves the $\delta_{\mathbf{k}\mathbf{k}'}\delta_{nn'}$ metric, i.e., the standard inner product — this is precisely the significance of the orthonormalization result proved in §3.3. If the metric were not the identity, the problem would become a generalized eigenvalue problem.)

Taking the variation of $R[\boldsymbol{\psi}]$ with respect to $\psi^{n*}_\mathbf{k}$ and setting $\delta R / \delta \psi^{n*}_\mathbf{k} = 0$, using the quotient rule:

$$
\frac{\delta}{\delta \psi^{n*}_\mathbf{k}} \left( \frac{\boldsymbol{\psi}^\dagger \mathcal{H}^{\text{eff}} \boldsymbol{\psi}}{\boldsymbol{\psi}^\dagger \boldsymbol{\psi}} \right) = 0
$$

$$
\Longrightarrow \quad \frac{1}{\boldsymbol{\psi}^\dagger \boldsymbol{\psi}} \left[ \sum_{\mathbf{p}, n'} \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}\, \psi^{n'}_\mathbf{p} - R[\boldsymbol{\psi}] \cdot \psi^n_\mathbf{k} \right] = 0
$$

Since $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} \neq 0$, the expression in brackets must vanish. Denoting the value of the Rayleigh quotient at the stationary point by $\varepsilon_\mu$, we obtain the standard **eigenvalue problem**:

$$
\boxed{\sum_{\mathbf{p}, n'} \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\, \psi^{n'}_\mathbf{p}(\mu, \mathbf{q}) = \varepsilon_\mu(\mathbf{q})\, \psi^n_\mathbf{k}(\mu, \mathbf{q})}
$$

This is the **Bethe-Salpeter equation** in the Tamm-Dancoff approximation. Its physical content is:

- Diagonalizing the matrix $\mathcal{H}^{\text{eff}}(\mathbf{q})$ yields eigenvalues $\varepsilon_\mu(\mathbf{q})$, the excitation energies of each collective mode $\mu$ at momentum $\mathbf{q}$
- The corresponding eigenvectors $\psi^n_\mathbf{k}(\mu, \mathbf{q})$ give the particle-hole pair wavefunction of that mode
- The **minimum** of the Rayleigh quotient gives the lowest-energy collective mode; **all stationary points** (saddle points) give the complete excitation spectrum

---

## 5. Complete Derivation of the Effective Hamiltonian $\mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})$

### 5.1 Decomposition Strategy

$$
\mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = (\mathcal{H}_0)^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + (\mathcal{H}_{\text{int}})^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
$$

We compute the contributions from the one-body part $[H_0, \cdot]$ and the two-body part $[H_{\text{int}}, \cdot]$ separately.

### 5.2 One-Body Part: $[H_0, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]$

Write $H_0$ in the quasiparticle basis:

$$
H_0 = \sum_{\mathbf{p}', mm'} \tilde{h}_{mm'}(\mathbf{p}')\, f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}'m'}
$$

where $\tilde{h}_{mm'}(\mathbf{p}') = \sum_{ab} U^*_{am}(\mathbf{p}')\, h^{ab}(\mathbf{p}')\, U_{bm'}(\mathbf{p}')$.

Using the commutation relation for fermionic bilinears (Appendix A):

$$
[f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}'m'},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}] = \delta_{\mathbf{p}', \mathbf{p}+\mathbf{q}}\,\delta_{m'n'}\, f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}, n_0} - \delta_{\mathbf{p}', \mathbf{p}}\,\delta_{m, n_0}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}'m'}
$$

Summing over $\mathbf{p}', m, m'$:

$$
[H_0,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}] = \sum_m \tilde{h}_{mn'}(\mathbf{p}+\mathbf{q})\, f^\dagger_{\mathbf{p}+\mathbf{q}, m}\, f_{\mathbf{p}, n_0} - \sum_{m'} \tilde{h}_{n_0 m'}(\mathbf{p})\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, m'}
$$

**Physical interpretation**: the first term scatters the particle ($n' \to m$); the second term scatters the hole ($n_0 \to m'$).

Taking the ground-state expectation value $\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n} \cdots |G\rangle$ and using the results of §3.3 ($n, n'$ unoccupied, $n_0$ occupied), the two terms respectively give:

$$
(\mathcal{H}_0)^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}} \left[ \tilde{h}_{nn'}(\mathbf{k}+\mathbf{q}) - \delta_{nn'}\, \tilde{h}_{n_0 n_0}(\mathbf{k}) \right]
$$

> **Note**: Here $\tilde{h}$ is the bare hopping Hamiltonian in the quasiparticle basis, which is generally **not diagonal**. The full mean-field Hamiltonian $\mathcal{H} = h + \Sigma$ is what diagonalizes. After combining the contributions from $H_0$ and $H_{\text{int}}$ (see §5.5), the mean-field part will fully diagonalize.

### 5.3 Two-Body Part: Expanding the Commutator

Write $H_{\text{int}}$ in the quasiparticle basis, transforming orbital indices to band indices with the $U$ matrix:

$$
H_{\text{int}} = \frac{1}{N}\sum_{\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}\sum_{m_1 m_2 m_3 m_4} \widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)\, f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}
$$

where $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2$, and the interaction in the quasiparticle basis is

$$
\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{p}_1)\, U_{bm_2}(\mathbf{p}_2)\, U^*_{cm_3}(\mathbf{p}_3)\, U_{dm_4}(\mathbf{p}_4)\; \widetilde{V}^{abcd}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)
$$

To compute $[H_{\text{int}},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]$, use

$$
[f^\dagger_1 f_2 f^\dagger_3 f_4,\; f^\dagger_5 f_6] = [f^\dagger_1 f_2,\; f^\dagger_5 f_6]\, f^\dagger_3 f_4 + f^\dagger_1 f_2\, [f^\dagger_3 f_4,\; f^\dagger_5 f_6]
$$

Each bilinear commutator produces two terms (Appendix A), yielding a total of **4 terms**, each containing four fermion operators:

**Part I**: $[f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]$

- Term (I-a): $\delta_{\mathbf{p}_2, \mathbf{p}+\mathbf{q}}\,\delta_{m_2,n'}$, producing $f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}, n_0}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}$
- Term (I-b): $\delta_{\mathbf{p}_1, \mathbf{p}}\,\delta_{m_1,n_0}$, producing $-f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}$

**Part II**: $[f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4},\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0}]$

- Term (II-a): $\delta_{\mathbf{p}_4, \mathbf{p}+\mathbf{q}}\,\delta_{m_4,n'}$, producing $f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}, n_0}$
- Term (II-b): $\delta_{\mathbf{p}_3, \mathbf{p}}\,\delta_{m_3,n_0}$, producing $-f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}_2,m_2}\, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}_4,m_4}$

### 5.4 Ground-State Expectation Values

For each of the four terms, we need to evaluate six-operator expectation values of the form

$$
\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [\text{4 } f/f^\dagger \text{ operators}]\; |G\rangle
$$

On the Slater-determinant ground state, Wick's theorem decomposes these into all possible full contractions (three pairs), with the contraction rule

$$
\langle G|\, f^\dagger_{\mathbf{p}m}\, f_{\mathbf{p}'m'}\, |G\rangle = \delta_{\mathbf{p}\mathbf{p}'}\,\delta_{mm'}\, \bar{n}_m(\mathbf{p})
$$

**Key selection rules**:

- The external $f_{\mathbf{k}+\mathbf{q}, n}$ ($n \in$ unocc): $\bar{n}_n = 0$, so it cannot pair with any $f^\dagger$ through a Wick contraction — it must "annihilate" against an $f^\dagger$ with matching momentum and band index via the anticommutation relation
- The external $f^\dagger_{\mathbf{k}, n_0}$ ($n_0 \in$ occ): $\bar{n}_{n_0} = 1$, so it can pair normally with a matching $f$

Therefore, in the six-operator expectation value, $f_{\mathbf{k}+\mathbf{q}, n}$ and $f^\dagger_{\mathbf{k}, n_0}$ must respectively match with the corresponding $f^\dagger$ and $f$ among the four internal operators via delta functions. The remaining internal pair then contracts with the Fermi sea.

**Analysis of the four terms** (Term (I-a) detailed; the others are analogous):

#### Term (I-a): Fix $\mathbf{p}_2 = \mathbf{p}+\mathbf{q}$, $m_2 = n'$

Six-operator structure:

$$
\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\, f^\dagger_{\mathbf{p}_1,m_1}\, f_{\mathbf{p}, n_0}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}\, |G\rangle
$$

Three classes of surviving contractions:

**Contraction (α)**: $f_{\mathbf{k}+\mathbf{q}, n} \leftrightarrow f^\dagger_{\mathbf{p}_1,m_1}$, $f^\dagger_{\mathbf{k}, n_0} \leftrightarrow f_{\mathbf{p}, n_0}$, $f^\dagger_{\mathbf{p}_3,m_3} \leftrightarrow f_{\mathbf{p}_4,m_4}$

Requires $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$; $\mathbf{k} = \mathbf{p}$; $\mathbf{p}_3 = \mathbf{p}_4$, $m_3 = m_4$, $\bar{n}_{m_3} \neq 0$.

The sum over $m_3 \in$ occ produces a density matrix — this is part of the Hartree-Fock self-energy correction to the particle energy.

**Contraction (β)**: $f_{\mathbf{k}+\mathbf{q}, n} \leftrightarrow f^\dagger_{\mathbf{p}_3,m_3}$, $f^\dagger_{\mathbf{k}, n_0} \leftrightarrow f_{\mathbf{p}, n_0}$, $f^\dagger_{\mathbf{p}_1,m_1} \leftrightarrow f_{\mathbf{p}_4,m_4}$

Requires $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $m_3 = n$; $\mathbf{k} = \mathbf{p}$; $\mathbf{p}_1 = \mathbf{p}_4$, $m_1 = m_4$, $\bar{n}_{m_1} \neq 0$.

This is also a self-energy correction type.

**Contraction (γ)**: $f_{\mathbf{k}+\mathbf{q}, n} \leftrightarrow f^\dagger_{\mathbf{p}_1,m_1}$, $f^\dagger_{\mathbf{k}, n_0} \leftrightarrow f_{\mathbf{p}_4,m_4}$, $f^\dagger_{\mathbf{p}_3,m_3} \leftrightarrow f_{\mathbf{p}, n_0}$

Requires $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$; $\mathbf{p}_4 = \mathbf{k}$, $m_4 = n_0$; $\mathbf{p}_3 = \mathbf{p}$, $m_3 = n_0$, $\bar{n}_{n_0} = 1$.

This is a **genuine particle-hole scattering** term — it does not involve a density-matrix summation over the Fermi sea, but rather describes the particle and hole scattering through the interaction.

After performing the analogous analysis for all four terms (I-a), (I-b), (II-a), (II-b), the contributions fall into two classes.

### 5.5 Combining the Results

**(A) Mean-field part**: From contractions of type (α) and (β) (one external operator pair matched + one internal pair contracted with the Fermi sea). These terms exactly reproduce the Hartree-Fock self-energy corrections to the particle and hole. Combined with the one-body part from §5.2, $\tilde{h}$ is replaced by the full mean-field Hamiltonian $\tilde{\mathcal{H}} = \tilde{h} + \tilde{\Sigma}$, which fully diagonalizes to the mean-field energies:

$$
(\mathcal{H}_{\text{MF}})^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right)
$$

This is the free particle-hole pair energy: particle energy minus hole energy.

**(B) Particle-hole interaction kernel**: From contractions of type (γ) (external particle operator paired with an internal $f^\dagger$, external hole operator paired with an internal $f$, no Fermi-sea contraction). These are the genuine residual interactions that cannot be absorbed into the mean field.

### 5.6 Two Topologies of the Particle-Hole Interaction Kernel

**(B1) Direct channel (density channel)** $\mathcal{K}^{\text{d}}$:

From the type-(γ) contractions in terms (I-a) and (II-a). Physical process: the interaction line connects the particle and hole lines, each scattering without exchange.

$$
\boxed{
\mathcal{K}^{\text{d},nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = -\frac{1}{N}\sum_{abcd} U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{bn'}(\mathbf{p}+\mathbf{q})\, U^*_{cn_0}(\mathbf{p})\, U_{dn_0}(\mathbf{k})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q}, \mathbf{p}+\mathbf{q}, \mathbf{p})
}
$$

Feynman diagram:

```
  particle: (p+q, n') ──[V]──▸ (k+q, n)
                            |
  hole:     (k, n₀)   ──[V]──▸ (p, n₀)
```

Momentum conservation check: in $\widetilde{V}^{abcd}(\mathbf{k}_1, \mathbf{k}_2, \mathbf{k}_3)$, $\mathbf{k}_1 = \mathbf{k}+\mathbf{q}$ (corresponding to $c^\dagger_a$, particle out), $\mathbf{k}_2 = \mathbf{p}+\mathbf{q}$ (corresponding to $c_b$, particle in), $\mathbf{k}_3 = \mathbf{p}$ (corresponding to $c^\dagger_c$, hole out), $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2 = \mathbf{k}$ (corresponding to $c_d$, hole in). ✓

**(B2) Exchange channel** $\mathcal{K}^{\text{x}}$:

From the line-crossing contractions in terms (I-a) and (II-b) (or (I-b) and (II-a)).

$$
\boxed{
\mathcal{K}^{\text{x},nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = +\frac{1}{N}\sum_{abcd} U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{dn_0}(\mathbf{k})\, U^*_{cn_0}(\mathbf{p})\, U_{bn'}(\mathbf{p}+\mathbf{q})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q}, \mathbf{k}, \mathbf{p})
}
$$

Feynman diagram (line crossing):

```
  particle: (p+q, n') ──╲    ╱──▸ (k+q, n)
                          ╲  ╱
                           ╳  [V]
                          ╱  ╲
  hole:     (k, n₀)   ──╱    ╲──▸ (p, n₀)
```

Momentum conservation check: $\mathbf{k}_1 = \mathbf{k}+\mathbf{q}$ (corresponding to $c^\dagger_a$, particle out), $\mathbf{k}_2 = \mathbf{k}$ (corresponding to $c_b$, hole in), $\mathbf{k}_3 = \mathbf{p}$ (corresponding to $c^\dagger_c$, hole out), $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2 = \mathbf{p}+\mathbf{q}$ (corresponding to $c_d$, particle in). ✓

> **Note on signs**: The direct channel carries a minus sign and the exchange channel carries a plus sign, arising from the relative fermion signs in the different contraction topologies. The specific signs depend on the operator ordering convention ($c^\dagger c c^\dagger c$) of the original Hamiltonian and the number of fermion line crossings in the Wick contractions.

### 5.7 Final Result

$$
\boxed{
\mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right) + \mathcal{K}^{\text{d},nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + \mathcal{K}^{\text{x},nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
}
$$

where:

- **First term**: free particle-hole pair energy (entirely from the mean field)
- **Second term**: direct channel kernel $\mathcal{K}^{\text{d}}$, with $\widetilde{V}$ three-momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$
- **Third term**: exchange channel kernel $\mathcal{K}^{\text{x}}$, with $\widetilde{V}$ three-momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$

In the general case both terms contribute. Depending on the symmetries of a specific model, one of them may vanish.

---

## 6. Relation to the Tamm-Dancoff and Random Phase Approximations

### 6.1 Tamm-Dancoff Approximation (TDA)

The eigenvalue problem derived above,

$$
\sum_{\mathbf{p}, n'} \mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})\, \psi^{n'}_\mathbf{p}(\mu, \mathbf{q}) = \varepsilon_\mu(\mathbf{q})\, \psi^n_\mathbf{k}(\mu, \mathbf{q})
$$

which includes only particle-hole ($ph$) excitations, is equivalent to the **Tamm-Dancoff approximation (TDA)**.

### 6.2 Full RPA

The full **random phase approximation (RPA)** additionally includes hole-particle ($hp$) backward excitations:

$$
|\mu, \mathbf{q}\rangle_{\text{RPA}} = \sum_{\mathbf{k}, n} \left[ X^n_\mathbf{k}\, f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0} - Y^n_\mathbf{k}\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n} \right] |G\rangle
$$

This leads to the standard matrix form of the RPA equation:

$$
\begin{pmatrix} A & B \\ -B^* & -A^* \end{pmatrix} \begin{pmatrix} X \\ Y \end{pmatrix} = \varepsilon \begin{pmatrix} X \\ Y \end{pmatrix}
$$

where the $A$ matrix is precisely $\mathcal{H}^{\text{eff}}$ (i.e., the TDA matrix), and the $B$ matrix describes ground-state correlations. The Tamm-Dancoff approximation corresponds to $B = 0$.

### 6.3 Collective Excitation Dispersion

The eigenvalues $\varepsilon_\mu(\mathbf{q})$ give the dispersion relation of collective excitations, with $\mu$ labeling different excitation bands. If the system spontaneously breaks a continuous symmetry, the Goldstone theorem guarantees the existence of a gapless collective mode at $\mathbf{q} \to 0$.

---

## 7. Summary: Complete Computational Workflow from the Original Hamiltonian to the Excitation Spectrum

1. **Hartree-Fock self-consistent calculation**:
   - Initialize the density matrix $G^{ab}(\mathbf{k})$
   - Construct $\mathcal{H}^{ab}(\mathbf{k}) = h^{ab}(\mathbf{k}) + \Sigma^{ab}(\mathbf{k})$
   - Diagonalize to obtain $E^n_\mathbf{k}$, $U(\mathbf{k})$
   - Update $G^{ab}(\mathbf{k})$, iterate until convergence

2. **Construct the effective Hamiltonian**:
   - For each $\mathbf{q}$, build the matrix $\mathcal{H}^{nn'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})$
   - Diagonal part: $E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}$
   - Off-diagonal part: $\mathcal{K}^{\text{d}} + \mathcal{K}^{\text{x}}$ (using the $U$ matrix to transform from orbital space to band space)

3. **Diagonalize the effective Hamiltonian**:
   - $\mathcal{H}^{\text{eff}}\, \psi = \varepsilon\, \psi$
   - Eigenvalues $\varepsilon_\mu(\mathbf{q})$ give the excitation spectrum
   - Eigenvectors $\psi^n_\mathbf{k}(\mu, \mathbf{q})$ give the excited-state wavefunctions

4. **Physical observables**:
   - Dynamic structure factor: $S(\mathbf{q}, \omega) = \sum_\mu |\langle \mu, \mathbf{q} | \hat{A}_\mathbf{q} | G \rangle|^2\, \delta(\omega - \varepsilon_\mu(\mathbf{q}))$
   - where $\hat{A}_\mathbf{q}$ is the operator that couples to the excitation of interest (e.g., spin operator, density operator); its matrix elements are computed via $U$ and $\psi$

---

## Appendix A: Useful Commutation Relations

For fermionic bilinear operators:

$$
[f^\dagger_\alpha f_\beta,\; f^\dagger_\gamma f_\delta] = \delta_{\beta\gamma}\, f^\dagger_\alpha f_\delta - \delta_{\alpha\delta}\, f^\dagger_\gamma f_\beta
$$

Derivation: using $\{f_\alpha, f^\dagger_\beta\} = \delta_{\alpha\beta}$,

$$
\begin{aligned}
f^\dagger_\alpha f_\beta f^\dagger_\gamma f_\delta &= f^\dagger_\alpha (\delta_{\beta\gamma} - f^\dagger_\gamma f_\beta) f_\delta = \delta_{\beta\gamma}\, f^\dagger_\alpha f_\delta - f^\dagger_\alpha f^\dagger_\gamma f_\beta f_\delta \\
f^\dagger_\gamma f_\delta f^\dagger_\alpha f_\beta &= \delta_{\delta\alpha}\, f^\dagger_\gamma f_\beta - f^\dagger_\gamma f^\dagger_\alpha f_\delta f_\beta = \delta_{\delta\alpha}\, f^\dagger_\gamma f_\beta - f^\dagger_\alpha f^\dagger_\gamma f_\beta f_\delta
\end{aligned}
$$

(The last step uses $f^\dagger_\gamma f^\dagger_\alpha = -f^\dagger_\alpha f^\dagger_\gamma$ and $f_\delta f_\beta = -f_\beta f_\delta$.) Subtracting the two expressions yields the result.

## Appendix B: Brillouin's Theorem

The matrix element of the full Hamiltonian $H$ between the Hartree-Fock ground state $|G\rangle$ and any single particle-hole excited state vanishes:

$$
\langle G | H | f^\dagger_{\mathbf{k}+\mathbf{q}, n}\, f_{\mathbf{k}, n_0} | G \rangle = 0 \quad (n \in \text{unocc}, \; n_0 \in \text{occ})
$$

This is a direct consequence of the Hartree-Fock equations (the self-consistency condition). Physically, the Hartree-Fock ground state is the "optimal" single Slater determinant, whose energy cannot be lowered by any single particle-hole excitation.

This theorem guarantees the validity of recasting $\varepsilon$ in commutator form in §4.2 ($H|G\rangle = E_G|G\rangle$ is effective within the single particle-hole subspace).

## Appendix C: $U$-Matrix Transformation and the Orbital-to-Band Conversion of $\widetilde{V}$

The interaction $\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)$ in the original orbital space is transformed to the band basis via the $U$ matrix:

$$
\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{k}_1)\, U_{bm_2}(\mathbf{k}_2)\, U^*_{cm_3}(\mathbf{k}_3)\, U_{dm_4}(\mathbf{k}_4)\; \widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)
$$

where $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2$.

When constructing $\mathcal{K}^{\text{d}}$ and $\mathcal{K}^{\text{x}}$, only the components of $\widetilde{V}$ with specific band indices ($n, n', n_0$) are needed. In practice, there is no need to fully transform the entire tensor — one can directly contract with the relevant columns of the $U$ matrix in orbital space.
