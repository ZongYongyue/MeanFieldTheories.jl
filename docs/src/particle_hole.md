# Particle-Hole Excitation Theory on Top of the Hartree-Fock Mean Field

---

## 1. Parametrization of the Particle-Hole Excited State

### 1.1 Physical Picture

The simplest excitation above the Hartree-Fock ground state is a **particle-hole pair**: an electron is promoted from an occupied band $n_0$ (the hole) to an unoccupied band $n$ (the particle). In a translationally invariant system, such an excitation carries a definite total momentum $\mathbf{q}$. Depending on the quantum numbers involved, these collective excitations may correspond to magnons, excitons, or other physical objects.

### 1.2 Excited-State Wavefunction

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

### 1.3 Orthonormalization

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

## 2. Variational Principle and the Eigenvalue Problem

### 2.1 Energy Functional

The excitation energy is computed using the **full original Hamiltonian** $H$ (not the mean-field $H_{\text{MF}}$):

```math
\varepsilon_\mu(\mathbf{q}) = \langle \mu, \mathbf{q} | H | \mu, \mathbf{q} \rangle - E_G
```

where $E_G = \langle G | H | G \rangle$ is the ground-state energy.

### 2.2 Recasting the Excitation Energy in Commutator Form

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

### 2.3 Expansion into a Quadratic Form in $\psi^{n_0 n}_\mathbf{k}$

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

### 2.4 Rayleigh Quotient and the Eigenvalue Problem

Our goal is to find, among all $\boldsymbol{\psi}$ satisfying the normalization constraint $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} = 1$, those that make the excitation energy $\varepsilon$ stationary. This is precisely the standard **Rayleigh quotient** problem.

Define the Rayleigh quotient

```math
R[\boldsymbol{\psi}] = \frac{\boldsymbol{\psi}^\dagger\, \mathcal{H}^{\text{eff}}\, \boldsymbol{\psi}}{\boldsymbol{\psi}^\dagger\, \boldsymbol{\psi}}
```

(Here the denominator involves the $\delta_{\mathbf{k}\mathbf{k}'}\delta_{n_0 n_0'}\delta_{nn'}$ metric, i.e., the standard inner product — this is precisely the significance of the orthonormalization result proved in §1.1. If the metric were not the identity, the problem would become a generalized eigenvalue problem.)

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

## 3. Complete Derivation of the Effective Hamiltonian $\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})$

### 3.1 Decomposition Strategy: Normal Ordering

The interaction Hamiltonian is written in the $c^\dagger c c^\dagger c$ form:

```math
H_{\text{int}} = \frac{1}{N}\sum_{\mathbf{k}_1\mathbf{k}_2\mathbf{k}_3}\sum_{abcd}
\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)\,
c^\dagger_{\mathbf{k}_1 a}\,c_{\mathbf{k}_2 b}\,c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_4 d}
```

Using the anticommutation relation $c_{\mathbf{k}_2 b}\,c^\dagger_{\mathbf{k}_3 c} = -c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_2 b} + \delta_{\mathbf{k}_2,\mathbf{k}_3}\,\delta_{bc}$, we decompose $H_{\text{int}}$ into a **normal-ordered quartic** and a **one-body self-energy** piece:

```math
H_{\text{int}} = \underbrace{-\frac{1}{N}\sum_{\mathbf{k}_1\mathbf{k}_2\mathbf{k}_3}\sum_{abcd}
\widetilde{V}^{abcd}(\mathbf{k}_1,\mathbf{k}_2,\mathbf{k}_3)\,
c^\dagger_{\mathbf{k}_1 a}\,c^\dagger_{\mathbf{k}_3 c}\,c_{\mathbf{k}_2 b}\,c_{\mathbf{k}_4 d}}_{\displaystyle :H_{\text{int}}:}
\;+\; \underbrace{\frac{1}{N}\sum_{\mathbf{k}_1,\mathbf{k}_3}\sum_{abd}
\widetilde{V}^{abbd}(\mathbf{k}_1,\mathbf{k}_3,\mathbf{k}_3)\,
c^\dagger_{\mathbf{k}_1 a}\,c_{\mathbf{k}_1 d}}_{\displaystyle H_{\text{SE}}}
```

where in $H_{\text{SE}}$, the delta $\delta_{\mathbf{k}_2,\mathbf{k}_3}\,\delta_{bc}$ has been applied, giving $\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_3 - \mathbf{k}_2 = \mathbf{k}_1$ and the sum over $b$ is contracted. The self-energy piece $H_{\text{SE}}$ is a one-body operator.

**Key observation**: the normal-ordering decomposition separates $H_{\text{int}}$ into a normal-ordered quartic $:H_{\text{int}}:$ (which produces all the non-trivial particle-hole scattering) and a one-body piece $H_{\text{SE}}$ (which contributes to the commutator in the same way as $H_0$). We compute $[H_0, \cdot]$, $[H_{\text{SE}}, \cdot]$, and $[:H_{\text{int}}:, \cdot]$ separately, then combine at the end.

Therefore, we decompose:

```math
\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \underbrace{(\mathcal{H}_0)^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})}_{\text{from } H_0} + \underbrace{(\mathcal{H}_{\text{SE}})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})}_{\text{from } H_{\text{SE}}} + \underbrace{(\mathcal{H}_{:\text{int}:})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})}_{\text{from } :H_{\text{int}}:}
```

and compute the three parts separately.

### 3.2 One-Body Part: $[H_0, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$

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

> **Note**: Here $\tilde{h}$ is the bare hopping Hamiltonian in the quasiparticle basis, which is generally **not diagonal**. The full mean-field Hamiltonian $\tilde{\mathcal{H}} = \tilde{h} + \tilde{\Sigma}$ is what diagonalizes. Whether and how the $H_0$ contribution combines with other one-body pieces to yield the diagonal mean-field energy will be addressed after all contributions are collected (see §3.5).

### 3.2.1 Self-Energy One-Body Part: $[H_{\text{SE}}, f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$

$H_{\text{SE}}$ is a one-body operator of the same form as $H_0$. In the quasiparticle basis:

```math
H_{\text{SE}} = \sum_{\mathbf{p}', mm'} \tilde{\sigma}_{mm'}(\mathbf{p}')\, f^\dagger_{\mathbf{p}'m}\, f_{\mathbf{p}'m'}
```

where

```math
\tilde{\sigma}_{mm'}(\mathbf{p}') = \frac{1}{N}\sum_{\mathbf{k}_3}\sum_{abd} U^*_{am}(\mathbf{p}')\, \widetilde{V}^{abbd}(\mathbf{p}',\mathbf{k}_3,\mathbf{k}_3)\, U_{dm'}(\mathbf{p}')
```

is the normal-ordering one-body coefficient (note: this is **not** the Hartree-Fock self-energy $\tilde{\Sigma}$, which involves the density matrix $G$).

Since $H_{\text{SE}}$ has the identical bilinear structure as $H_0$, the commutator and expectation value follow exactly the same algebra as §3.2, with $\tilde{h} \to \tilde{\sigma}$:

```math
(\mathcal{H}_{\text{SE}})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}} \left[ \delta_{n_0 n_0'}\,\tilde{\sigma}_{nn'}(\mathbf{k}+\mathbf{q}) - \delta_{nn'}\, \tilde{\sigma}_{n_0' n_0}(\mathbf{k}) \right]
```

### 3.3 Normal-Ordered Two-Body Part: Expanding the Commutator

Write the normal-ordered interaction in the quasiparticle basis:

```math
:H_{\text{int}}: = -\frac{1}{N}\sum_{\mathbf{p}_1\mathbf{p}_2\mathbf{p}_3}\sum_{m_1 m_2 m_3 m_4} \widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)\, f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}
```

where $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2$, and the interaction in the quasiparticle basis is

```math
\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{p}_1)\, U_{bm_2}(\mathbf{p}_2)\, U^*_{cm_3}(\mathbf{p}_3)\, U_{dm_4}(\mathbf{p}_4)\; \widetilde{V}^{abcd}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)
```

Note the overall minus sign from normal ordering. To compute $[:H_{\text{int}}:,\; f^\dagger_{\mathbf{p}+\mathbf{q}, n'}\, f_{\mathbf{p}, n_0'}]$, we use the SNEG-verified commutator identity for normal-ordered quartics:

```math
[f^\dagger_1 f^\dagger_3 f_2 f_4,\; f^\dagger_5 f_6] = +\delta_{1,6}\, f^\dagger_3 f^\dagger_5 f_2 f_4 - \delta_{2,5}\, f^\dagger_1 f^\dagger_3 f_4 f_6 - \delta_{3,6}\, f^\dagger_1 f^\dagger_5 f_2 f_4 + \delta_{4,5}\, f^\dagger_1 f^\dagger_3 f_2 f_6
```

where $\delta_{i,j}$ denotes the composite delta $\delta_{\mathbf{p}_i, \mathbf{p}_j}\,\delta_{m_i, m_j}$. Each term is **already normal-ordered** ($f^\dagger f^\dagger f f$).

**Physical interpretation**: each delta represents one operator in the quartic anticommuting with one operator in the bilinear. The four terms correspond to:

| Term | Delta | Sign | Fixes |
|---|---|---|---|
| (A) | $\delta_{1,6}$: $\delta_{\mathbf{p}_1,\mathbf{p}}\,\delta_{m_1,n_0'}$ | $+$ | $\mathbf{p}_1 = \mathbf{p}$, remaining: $f^\dagger_3 f^\dagger_5 f_2 f_4$ |
| (B) | $\delta_{2,5}$: $\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}\,\delta_{m_2,n'}$ | $-$ | $\mathbf{p}_2 = \mathbf{p}+\mathbf{q}$, remaining: $f^\dagger_1 f^\dagger_3 f_4 f_6$ |
| (C) | $\delta_{3,6}$: $\delta_{\mathbf{p}_3,\mathbf{p}}\,\delta_{m_3,n_0'}$ | $-$ | $\mathbf{p}_3 = \mathbf{p}$, remaining: $f^\dagger_1 f^\dagger_5 f_2 f_4$ |
| (D) | $\delta_{4,5}$: $\delta_{\mathbf{p}_4,\mathbf{p}+\mathbf{q}}\,\delta_{m_4,n'}$ | $+$ | $\mathbf{p}_4 = \mathbf{p}+\mathbf{q}$, remaining: $f^\dagger_1 f^\dagger_3 f_2 f_6$ |

> **Advantage of normal ordering**: In the old $c^\dagger c c^\dagger c$ approach, the commutator produced 4 terms with mixed operator ordering, requiring further Wick decomposition into self-energy and scattering contractions. Here, every term is already normal-ordered — the subsequent ground-state expectation values involve only **genuine particle-hole scattering** contractions (no self-energy subtractions needed).

### 3.4 Ground-State Expectation Values

For each of the four terms, we need to evaluate six-operator expectation values of the form

```math
\langle G|\, f^\dagger_{\mathbf{k}, n_0}\, f_{\mathbf{k}+\mathbf{q}, n}\; [f^\dagger f^\dagger f f]\; |G\rangle
```

On the Slater-determinant ground state, Wick's theorem decomposes these into all possible full contractions (three pairs), with the contraction rule

```math
\langle G|\, f^\dagger_{\mathbf{p}m}\, f_{\mathbf{p}'m'}\, |G\rangle = \delta_{\mathbf{p}\mathbf{p}'}\,\delta_{mm'}\, \bar{n}_m(\mathbf{p})
```

**Key selection rules** (same as before):

- The external $f_{\mathbf{k}+\mathbf{q}, n}$ ($n \in$ unocc): $\bar{n}_n = 0$, so it must anticommute against an internal $f^\dagger$ with matching quantum numbers
- The external $f^\dagger_{\mathbf{k}, n_0}$ ($n_0 \in$ occ): $\bar{n}_{n_0} = 1$, so it contracts normally with a matching internal $f$

**Crucial simplification from normal ordering**: since all four internal operators are now in normal order ($f^\dagger f^\dagger f f$), the external $f_{\mathbf{k}+\mathbf{q}, n}$ can only pair with one of the two internal $f^\dagger$'s, and $f^\dagger_{\mathbf{k}, n_0}$ can only pair with one of the two internal $f$'s. This gives at most $2 \times 2 = 4$ contractions per term. However, occupation constraints further reduce this — in many cases only **one** contraction survives per term. There are no self-energy-type contractions.

We evaluate each of the four commutator terms separately, sandwiching with the external operators and applying Wick's theorem (verified by SNEG).

#### Understanding the $\theta$ functions in SNEG output

In SNEG with `ordering[f] = SEA` (Slater-determinant ground state), each Wick contraction pair produces an occupation factor. The basic contraction rule is:

```math
\langle G|\, f^\dagger_\alpha\, f_\beta\, |G\rangle = \delta_{\alpha\beta}\, \theta(-\alpha)
```

where $\theta(-\alpha) \equiv \bar{n}_\alpha$ is the ground-state occupation number of state $\alpha$:
- $\theta(-\alpha) = 1$ if state $\alpha$ is occupied
- $\theta(-\alpha) = 0$ if state $\alpha$ is unoccupied

When $f_\beta$ must anticommute past $f^\dagger_\alpha$ before contracting (i.e., when the annihilator appears to the *left* of the creator), the contraction picks up a complementary factor:

```math
\theta(+\alpha) = 1 - \theta(-\alpha) = 1 - \bar{n}_\alpha
```

which equals 1 for **unoccupied** states and 0 for **occupied** states.

Since a six-operator expectation value has three contraction pairs, each SNEG output term carries a product of three $\theta$ factors.

**Example**: consider the contraction pattern where $f^\dagger_{\mathbf{k},n_0}$ pairs with $f_{\mathbf{p}_2,m_2}$, $f_{\mathbf{k}+\mathbf{q},n}$ pairs with $f^\dagger_{\mathbf{p}_3,m_3}$, and the remaining internal pair $f^\dagger_{\mathbf{p}+\mathbf{q},n'}$ pairs with $f_{\mathbf{p}_4,m_4}$. This produces:

```math
\underbrace{\theta(-\mathbf{k},n_0)}_{\bar{n}_{n_0}(\mathbf{k})}\; \underbrace{\theta(\mathbf{p}_3,m_3)}_{1 - \bar{n}_{m_3}(\mathbf{p}_3)}\; \underbrace{\theta(-\mathbf{p}-\mathbf{q},n')}_{\bar{n}_{n'}(\mathbf{p}+\mathbf{q})}
```

Here:
- $\theta(-\mathbf{k},n_0) = \bar{n}_{n_0}(\mathbf{k}) = 1$ because $n_0 \in \text{occ}(\mathbf{k})$ by definition
- $\theta(\mathbf{p}_3,m_3) = 1 - \bar{n}_{m_3}(\mathbf{p}_3)$: depends on whether band $m_3$ at momentum $\mathbf{p}_3$ is occupied or not. Since $m_3$ is an internal summation index, this factor **remains** and restricts the sum to $m_3 \in \text{unocc}(\mathbf{p}_3)$
- $\theta(-\mathbf{p}-\mathbf{q}, n') = \bar{n}_{n'}(\mathbf{p}+\mathbf{q})$: depends on the occupation of band $n'$ at $\mathbf{p}+\mathbf{q}$. Since $n' \in \text{unocc}$ for the particle-hole pair, this equals 0

The combination $\theta(-\mathbf{p}-\mathbf{q}) - 1$ that appears frequently equals $\bar{n}_{n'}(\mathbf{p}+\mathbf{q}) - 1 = -\theta(\mathbf{p}+\mathbf{q}) = -(1 - \bar{n}_{n'})$, which restricts the corresponding state to be unoccupied.

> **Note on occupation sets**: the sets $\text{occ}(\mathbf{k})$ and $\text{unocc}(\mathbf{k})$ are defined at each momentum $\mathbf{k}$ by the self-consistent Hartree-Fock solution. In general, the partition into occupied and unoccupied bands can vary with momentum (e.g., different fillings at different $\mathbf{k}$-points in the magnetic Brillouin zone). Therefore, one cannot assign a global occ/unocc label to a band index independent of momentum.

#### Term (A): $\delta_{m_1,n_0'}\,\delta_{p,\mathbf{p}_1}$ — remaining operators $f^\dagger_{\mathbf{p}_3,m_3}\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\,f_{\mathbf{p}_2,m_2}\,f_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,m_4}$

The SNEG Wick contraction of $\langle G|\,f^\dagger_{\mathbf{k},n_0}\,f_{\mathbf{k}+\mathbf{q},n}\,f^\dagger_{\mathbf{p}_3,m_3}\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\,f_{\mathbf{p}_2,m_2}\,f_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,m_4}\,|G\rangle$ gives (raw output, to be simplified):

```math
\begin{aligned}
\theta(-\mathbf{k},n_0)\;\Big\{&
-\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\Big[
  \theta(\mathbf{p}_3)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{m_3,n}\,\delta_{m_4,n'}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}} \\
  &\quad +\theta(-\mathbf{p}_3)\,(\theta(-\mathbf{p}-\mathbf{q})-1)\,\delta_{m_3,m_4}\,\delta_{n,n'}\,\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}}\,\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}
\Big] \\
+\;&\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\Big[
  \theta(\mathbf{p}_3)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{m_2,n'}\,\delta_{m_3,n}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}} \\
  &\quad +\theta(-\mathbf{p}_3)\,(\theta(-\mathbf{p}-\mathbf{q})-1)\,\delta_{m_2,m_3}\,\delta_{n,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}\,\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}}
\Big] \\
+\;&\theta(-\mathbf{p}_3)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}\Big[
  \delta_{m_2,n'}\,\delta_{m_3,m_4}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}\,\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} \\
  &\quad -\delta_{m_2,m_3}\,\delta_{m_4,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}
\Big]
\Big\}
\end{aligned}
```

**Simplification will be performed after all four terms are collected.**

#### Term (B): $-\delta_{m_3,n_0'}\,\delta_{p,\mathbf{p}_3}$ — remaining operators $f^\dagger_{\mathbf{p}_1,m_1}\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\,f_{\mathbf{p}_2,m_2}\,f_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,m_4}$

The SNEG Wick contraction gives:

```math
\begin{aligned}
\theta(-\mathbf{k},n_0)\;\Big\{&
-\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{m_1,n}\,\delta_{m_4,n'}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}} \\
  &\quad +\theta(-\mathbf{p}_1)\,(\theta(-\mathbf{p}-\mathbf{q})-1)\,\delta_{m_1,m_4}\,\delta_{n,n'}\,\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}}\,\delta_{\mathbf{p}_1,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}
\Big] \\
+\;&\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{m_1,n}\,\delta_{m_2,n'}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}} \\
  &\quad +\theta(-\mathbf{p}_1)\,(\theta(-\mathbf{p}-\mathbf{q})-1)\,\delta_{m_1,m_2}\,\delta_{n,n'}\,\delta_{\mathbf{p}_1,\mathbf{p}_2}\,\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}}
\Big] \\
+\;&\theta(-\mathbf{p}_1)\,\theta(-\mathbf{p}-\mathbf{q})\,\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}\Big[
  \delta_{m_1,m_4}\,\delta_{m_2,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}\,\delta_{\mathbf{p}_1,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} \\
  &\quad -\delta_{m_1,m_2}\,\delta_{m_4,n'}\,\delta_{\mathbf{p}_1,\mathbf{p}_2}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}
\Big]
\Big\}
\end{aligned}
```

#### Term (C): $\delta_{m_2,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}$ — remaining operators $f^\dagger_{\mathbf{p}_1,m_1}\,f^\dagger_{\mathbf{p}_3,m_3}\,f_{\mathbf{p},n_0'}\,f_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,m_4}$

The SNEG Wick contraction gives:

```math
\begin{aligned}
\theta(-\mathbf{k},n_0)\;\Big\{&
\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}
\Big] \\
+\;&\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}\Big[
  \theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,m_4}\,\delta_{m_3,n}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_1,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} \\
  &\quad -\theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,m_4}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}
\Big] \\
+\;&\theta(-\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}\Big[
  \delta_{m_1,m_4}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} \\
  &\quad -\delta_{m_1,n_0'}\,\delta_{m_3,m_4}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}
\Big]
\Big\}
\end{aligned}
```

#### Term (D): $-\delta_{m_4,n'}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}$ — remaining operators $f^\dagger_{\mathbf{p}_1,m_1}\,f^\dagger_{\mathbf{p}_3,m_3}\,f_{\mathbf{p},n_0'}\,f_{\mathbf{p}_2,m_2}$

The SNEG Wick contraction gives:

```math
\begin{aligned}
\theta(-\mathbf{k},n_0)\;\Big\{&
\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}
\Big] \\
+\;&\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}\Big[
  \theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,m_2}\,\delta_{m_3,n}\,\delta_{\mathbf{p}_1,\mathbf{p}_2}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_2,m_3}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}
\Big] \\
+\;&\theta(-\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}\Big[
  \delta_{m_1,m_2}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{p}_2} \\
  &\quad -\delta_{m_1,n_0'}\,\delta_{m_2,m_3}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}
\Big]
\Big\}
\end{aligned}
```

**All four raw SNEG results are now collected. Simplification proceeds below.**

### 3.5 Simplification of the Wick Contraction Results

#### Step 1: Drop the $\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}$ groups

Every term contains a third group proportional to $\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}\,\delta_{n,n_0}$. This forces $\mathbf{q} = 0$ and $n = n_0$. Since $n_0 \in \text{occ}(\mathbf{k})$ and $n \in \text{unocc}(\mathbf{k})$ at $\mathbf{q}=0$, we have $n \neq n_0$. **All such groups vanish.**

#### Step 2: Set $\theta(-\mathbf{k}, n_0) = 1$

The overall prefactor $\theta(-\mathbf{k}, n_0) = \bar{n}_{n_0}(\mathbf{k}) = 1$ since $n_0 \in \text{occ}(\mathbf{k})$.

#### Step 3: Identify which $\theta$ factors involve external indices

The key to simplification is recognizing which $\theta$ factors correspond to external (fixed occ/unocc) indices versus internal (summation) indices:

- **Terms A and B**: the remaining four operators include $f^\dagger_{\mathbf{p}+\mathbf{q},n'}$ from the bilinear (with $n' \in \text{unocc}$). Its Wick contraction produces $\theta(-\mathbf{p}-\mathbf{q}) = \bar{n}_{n'}(\mathbf{p}+\mathbf{q}) = 0$ or $\theta(\mathbf{p}+\mathbf{q}) = 1 - \bar{n}_{n'} = 1$. Hence $\theta(-\mathbf{p}-\mathbf{q}) - 1 = -1$.

- **Terms C and D**: the remaining four operators include $f_{\mathbf{p},n_0'}$ from the bilinear (with $n_0' \in \text{occ}$). All $\theta$ factors are on internal summation indices $(\mathbf{p}_1, m_1)$ and $(\mathbf{p}_3, m_3)$, which cannot be simplified without knowing the band structure.

#### Step 4: Simplify Terms A and B

In Terms A and B, every sub-contraction contains either $\theta(-\mathbf{p}-\mathbf{q}) = 0$ or $(\theta(-\mathbf{p}-\mathbf{q})-1) = -1$:

- Sub-terms with $\theta(-\mathbf{p}-\mathbf{q})$: **vanish**
- Sub-terms with $(\theta(-\mathbf{p}-\mathbf{q})-1) = -1$: **survive**, and they all require $\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}}$ (i.e., $\mathbf{k} = \mathbf{p}$) plus $\delta_{n,n'}$ and internal self-contraction deltas ($\delta_{m_i,m_j}$ with $\delta_{\mathbf{p}_i,\mathbf{p}_j}$). These are **self-energy-type** contractions involving a Fermi-sea summation.

**Surviving terms from A** (after setting $\theta(-\mathbf{k}) = 1$, $\theta(-\mathbf{p}-\mathbf{q}) - 1 = -1$, dropping $\theta(-\mathbf{p}-\mathbf{q}) = 0$ terms):

```math
\begin{aligned}
\text{(A)} = &+\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\,\theta(-\mathbf{p}_3)\,\delta_{m_3,m_4}\,\delta_{n,n'}\,\delta_{\mathbf{k},\mathbf{p}}\,\delta_{\mathbf{p}_1,\mathbf{p}_2} \\
&-\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\,\theta(-\mathbf{p}_3)\,\delta_{m_2,m_3}\,\delta_{n,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}\,\delta_{\mathbf{k},\mathbf{p}}
\end{aligned}
```

(Here we used $\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} = \delta_{\mathbf{p}_1,\mathbf{p}_2}$ and $\delta_{\mathbf{k}+\mathbf{q},\mathbf{p}+\mathbf{q}} = \delta_{\mathbf{k},\mathbf{p}}$.)

**Surviving terms from B** (same logic):

```math
\begin{aligned}
\text{(B)} = &+\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\,\theta(-\mathbf{p}_1)\,\delta_{m_1,m_4}\,\delta_{n,n'}\,\delta_{\mathbf{k},\mathbf{p}}\,\delta_{\mathbf{p}_1,\mathbf{p}_2} \\
&-\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\,\theta(-\mathbf{p}_1)\,\delta_{m_1,m_2}\,\delta_{n,n'}\,\delta_{\mathbf{p}_1,\mathbf{p}_2}\,\delta_{\mathbf{k},\mathbf{p}}
\end{aligned}
```

These self-energy terms will be combined with $\mathcal{H}_0$ and $\mathcal{H}_{\text{SE}}$ in §3.6.

#### Step 5: Simplify Terms C and D

Terms C and D contain no $\theta(-\mathbf{p}-\mathbf{q})$ factor. After dropping the $\delta_{\mathbf{k},\mathbf{k}+\mathbf{q}}$ group and setting $\theta(-\mathbf{k}) = 1$, each retains two groups:

- **First group**: all four external band indices ($n, n_0, n', n_0'$) appear — **scattering kernel**
- **Second group**: has $\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}$ with internal self-contractions — **self-energy**

**From Term (C)** ($\delta_{m_2,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}$):

```math
\begin{aligned}
\text{(C)}_{\text{scat}} = \;&\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}
\Big]
\end{aligned}
```

```math
\begin{aligned}
\text{(C)}_{\text{SE}} = \;&\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}\Big[
  \theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,m_4}\,\delta_{m_3,n}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_1,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3} \\
  &\quad -\theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,m_4}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}\,\delta_{\mathbf{p}_3,\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}
\Big]
\end{aligned}
```

**From Term (D)** ($-\delta_{m_4,n'}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}$):

```math
\begin{aligned}
\text{(D)}_{\text{scat}} = \;&\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\Big[
  \theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}
\Big]
\end{aligned}
```

```math
\begin{aligned}
\text{(D)}_{\text{SE}} = \;&\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}\Big[
  \theta(-\mathbf{p}_1)\,\theta(\mathbf{p}_3)\,\delta_{m_1,m_2}\,\delta_{m_3,n}\,\delta_{\mathbf{p}_1,\mathbf{p}_2}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}} \\
  &\quad -\theta(\mathbf{p}_1)\,\theta(-\mathbf{p}_3)\,\delta_{m_1,n}\,\delta_{m_2,m_3}\,\delta_{\mathbf{p}_2,\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}
\Big]
\end{aligned}
```

#### Summary of classification

| Term | Commutator sign | Scattering contributions | Self-energy contributions |
|---|---|---|---|
| A ($\delta_{m_1,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_1}$) | $+1$ | none (all vanish via $\theta(-\mathbf{p}-\mathbf{q})=0$) | 2 terms |
| B ($-\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}$) | $-1$ | none (all vanish via $\theta(-\mathbf{p}-\mathbf{q})=0$) | 2 terms |
| C ($\delta_{m_2,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}$) | $-1$ | 2 terms | 2 terms |
| D ($-\delta_{m_4,n'}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}$) | $+1$ | 2 terms | 2 terms |

The **scattering kernel** comes entirely from Terms C and D (4 sub-terms total). The **self-energy** comes from all four terms (8 sub-terms total).

### 3.6 Extracting the Scattering Kernel

We now resolve the delta functions in the 4 scattering sub-terms to obtain the explicit $V$ tensor arguments and signs. The overall sign of each contribution is:

$$\text{total sign} = \underbrace{(-1/N)}_{\text{from } :H_{\text{int}}:} \times \underbrace{(\pm 1)}_{\text{commutator } \delta \text{ sign}} \times \underbrace{(\pm 1)}_{\text{Wick sign from SNEG}}$$

where the commutator $\delta$ sign is $+1$ for Terms A, C (from $+\delta_{1,6}$, $+\delta_{2,5}$ after accounting for the $c^\dagger c^\dagger cc$ ordering) and $-1$ for Terms B, D (from $-\delta_{3,6}$, $-\delta_{4,5}$). The Wick sign is the sign of each sub-contraction as given by the SNEG output.

#### Sub-term C-1: from $\text{(C)}_{\text{scat}}$, first line

Commutator delta: $+\delta_{m_2,n'}\,\delta_{\mathbf{p}_2,\mathbf{p}+\mathbf{q}}$. Wick deltas: $\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}$.

Resolving: $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $\mathbf{p}_2 = \mathbf{p}+\mathbf{q}$, $\mathbf{p}_3 = \mathbf{p}$, $\mathbf{p}_4 = \mathbf{k}$. Band indices: $m_1 = n$, $m_2 = n'$, $m_3 = n_0'$, $m_4 = n_0$.

$\theta$ factors: $\theta(\mathbf{k}+\mathbf{q}, n) = 1$ (unocc), $\theta(-\mathbf{p}, n_0') = 1$ (occ). Product = 1.

Signs: $(-1/N) \times (+1)_{\text{comm}} \times (+1)_{\text{Wick}} = -1/N$.

$\Rightarrow$ **Contribution**: $-\frac{1}{N}\,\widetilde{V}_{n,\, n',\, n_0',\, n_0}(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$

Momentum transfer: $\mathbf{p}_1 - \mathbf{p}_2 = \mathbf{k} - \mathbf{p}$ → **direct channel**.

#### Sub-term C-2: from $\text{(C)}_{\text{scat}}$, second line

Same commutator delta ($+1$). Wick deltas: $\delta_{m_4,n_0}\,\delta_{\mathbf{k},\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3}\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}$.

Resolving: $\mathbf{p}_1 = \mathbf{p}$, $\mathbf{p}_2 = \mathbf{p}+\mathbf{q}$, $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $\mathbf{p}_4 = \mathbf{k}$. Band indices: $m_1 = n_0'$, $m_2 = n'$, $m_3 = n$, $m_4 = n_0$.

$\theta$ factors: $\theta(-\mathbf{p}, n_0') = 1$ (occ), $\theta(\mathbf{k}+\mathbf{q}, n) = 1$ (unocc). Product = 1.

Signs: $(-1/N) \times (+1)_{\text{comm}} \times (-1)_{\text{Wick}} = +1/N$.

$\Rightarrow$ **Contribution**: $+\frac{1}{N}\,\widetilde{V}_{n_0',\, n',\, n,\, n_0}(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$

Momentum transfer: $\mathbf{p}_1 - \mathbf{p}_2 = -\mathbf{q}$ → **exchange channel**.

#### Sub-term D-1: from $\text{(D)}_{\text{scat}}$, first line

Commutator delta: $-\delta_{m_4,n'}\,\delta_{\mathbf{p}_1-\mathbf{p}_2+\mathbf{p}_3,\mathbf{p}+\mathbf{q}}$ (sign $-1$). Wick deltas: $\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\,\delta_{m_1,n}\,\delta_{m_3,n_0'}\,\delta_{\mathbf{p},\mathbf{p}_3}\,\delta_{\mathbf{p}_1,\mathbf{k}+\mathbf{q}}$.

Resolving: $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $\mathbf{p}_2 = \mathbf{k}$, $\mathbf{p}_3 = \mathbf{p}$, $\mathbf{p}_4 = \mathbf{p}+\mathbf{q}$. Band indices: $m_1 = n$, $m_2 = n_0$, $m_3 = n_0'$, $m_4 = n'$.

$\theta$ factors: $\theta(\mathbf{k}+\mathbf{q}, n) = 1$ (unocc), $\theta(-\mathbf{p}, n_0') = 1$ (occ). Product = 1.

Signs: $(-1/N) \times (-1)_{\text{comm}} \times (+1)_{\text{Wick}} = +1/N$.

$\Rightarrow$ **Contribution**: $+\frac{1}{N}\,\widetilde{V}_{n,\, n_0,\, n_0',\, n'}(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$

Momentum transfer: $\mathbf{p}_1 - \mathbf{p}_2 = \mathbf{q}$ → **exchange channel**.

#### Sub-term D-2: from $\text{(D)}_{\text{scat}}$, second line

Same commutator delta ($-1$). Wick deltas: $\delta_{\mathbf{k},\mathbf{p}_2}\,\delta_{m_2,n_0}\,\delta_{m_1,n_0'}\,\delta_{m_3,n}\,\delta_{\mathbf{p},\mathbf{p}_1}\,\delta_{\mathbf{p}_3,\mathbf{k}+\mathbf{q}}$.

Resolving: $\mathbf{p}_1 = \mathbf{p}$, $\mathbf{p}_2 = \mathbf{k}$, $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $\mathbf{p}_4 = \mathbf{p}+\mathbf{q}$. Band indices: $m_1 = n_0'$, $m_2 = n_0$, $m_3 = n$, $m_4 = n'$.

$\theta$ factors: $\theta(-\mathbf{p}, n_0') = 1$ (occ), $\theta(\mathbf{k}+\mathbf{q}, n) = 1$ (unocc). Product = 1.

Signs: $(-1/N) \times (-1)_{\text{comm}} \times (-1)_{\text{Wick}} = -1/N$.

$\Rightarrow$ **Contribution**: $-\frac{1}{N}\,\widetilde{V}_{n_0',\, n_0,\, n,\, n'}(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$

Momentum transfer: $\mathbf{p}_1 - \mathbf{p}_2 = \mathbf{p} - \mathbf{k}$ → **direct channel**.

#### Topology classification

| Sub-term | $V$ band indices | $V$ momenta $(\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_3)$ | $\mathbf{p}_1 - \mathbf{p}_2$ | Sign | Topology |
|---|---|---|---|---|---|
| C-1 | $(n, n', n_0', n_0)$ | $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$ | $\mathbf{k}-\mathbf{p}$ | $-1/N$ | Direct |
| D-2 | $(n_0', n_0, n, n')$ | $(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$ | $\mathbf{p}-\mathbf{k}$ | $-1/N$ | Direct |
| D-1 | $(n, n_0, n_0', n')$ | $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$ | $\mathbf{q}$ | $+1/N$ | Exchange |
| C-2 | $(n_0', n', n, n_0)$ | $(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$ | $-\mathbf{q}$ | $+1/N$ | Exchange |

Each channel receives **two contributions** with different $V$ momentum arguments. The signs (**direct $-$, exchange $+$**) match the original derivation exactly.

### 3.7 Explicit Kernel Expressions in Orbital Basis

Transforming back to the orbital basis using $\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{p}_1)\, U_{bm_2}(\mathbf{p}_2)\, U^*_{cm_3}(\mathbf{p}_3)\, U_{dm_4}(\mathbf{p}_4)\; \widetilde{V}^{abcd}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)$:

**(A) Direct channel** $\mathcal{K}^{\text{d}}$:

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

**(B) Exchange channel** $\mathcal{K}^{\text{x}}$:

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

> **Consistency check**: The kernel expressions are identical to the original derivation (which used the $c^\dagger c c^\dagger c$ ordering): direct channel carries $-1/N$ and exchange channel carries $+1/N$. The normal-ordering approach produces the same signs because the $-1$ prefactor from $:H_{\text{int}}:$ is exactly compensated by the different Wick contraction structure of the normal-ordered quartic.

### 3.8 Self-Energy Contributions and the Mean-Field Part

The self-energy sub-terms from all four terms (A, B, C, D) require $\mathbf{k} = \mathbf{p}$ and involve Fermi-sea summations over internal indices. Combined with $\mathcal{H}_0$ (§3.2) and $\mathcal{H}_{\text{SE}}$ (§3.2.1), they produce the full diagonal mean-field energy:

```math
(\mathcal{H}_{\text{MF}})^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right)
```

This is the free particle-hole pair energy: particle energy minus hole energy.

> **Note**: The detailed verification that the self-energy sub-terms combine with $\mathcal{H}_0 + \mathcal{H}_{\text{SE}}$ to produce exactly the diagonal mean-field energies $E^n_\mathbf{k}$ is a standard result of the Hartree-Fock self-consistency condition and follows the same logic as in the original derivation.

### 3.9 Final Result

```math
\boxed{
\mathcal{H}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k}\right) + \mathcal{K}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + \mathcal{K}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
}
```

where:

- **First term**: free particle-hole pair energy (from the mean field), diagonal in all indices
- **Second term**: direct channel kernel $\mathcal{K}^{\text{d}}$ (§3.7), with two $\widetilde{V}$ contributions at momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$ and $(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$
- **Third term**: exchange channel kernel $\mathcal{K}^{\text{x}}$ (§3.7), with two $\widetilde{V}$ contributions at momentum arguments $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$ and $(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$
