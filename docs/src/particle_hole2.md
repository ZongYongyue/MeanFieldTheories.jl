# Particle-Hole and Hole-Particle Excitation Theory on Top of the Hartree-Fock Mean Field

---

## 1. Full Random Phase Approximation (RPA)

In the Tamm-Dancoff Approximation (TDA) discussed in the [particle-hole excitation theory](particle_hole.md), we assume the ground state is the absolute Hartree-Fock vacuum, and collective excitations consist purely of creating particle-hole pairs. However, in systems with strong correlations or spontaneous symmetry breaking (such as an antiferromagnetic state), the true ground state contains zero-point quantum fluctuations and is pre-mixed with "virtual" particle-hole pairs.

The full **Random Phase Approximation (RPA)** allows collective excitations not only to *create* particle-hole pairs (the forward process) but also to *annihilate* existing particle-hole pairs from the correlated ground state (the backward process).

### 1.1 The RPA Excitation Operator

In the TDA, the excitation operator $\hat{O}^\dagger$ contained only a single (forward) term. For RPA, we generalize it to include a backward (de-excitation) component:

```math
\hat{O}^\dagger_{\mu q} = \sum_{k, n_0, n} X^{n_0 n}_k(\mu, q)\, f^\dagger_{k+q, n}\, f_{k, n_0} - \sum_{k, n_0, n} Y^{n_0 n}_k(\mu, q)\, f^\dagger_{k, n_0}\, f_{k-q, n}
```

where $n_0 \in$ occ, $n \in$ unocc throughout:
- **Forward (excitation) part**: $f^\dagger_{k+q, n}\, f_{k, n_0}$ creates a particle (unocc $n$) at $k+q$ and a hole (occ $n_0$) at $k$. Net momentum: $+q$.
- **Backward (de-excitation) part**: $f^\dagger_{k, n_0}\, f_{k-q, n}$ creates an occupied-band particle at $k$ and removes an unoccupied-band particle at $k-q$. Net momentum: $+q$. This is **not** the Hermitian conjugate of the forward bilinear — rather, it carries the same net momentum $+q$ by shifting the momentum arguments while keeping the occ/unocc structure.

> **Why this momentum structure**: The backward bilinear must carry net momentum $+q$ to couple with the forward part through the momentum-conserving Hamiltonian. The choice $f^\dagger_{k, n_0}\, f_{k-q, n}$ achieves this: momentum in = $k-q$, momentum out = $k$, net = $+q$.

The Hermitian conjugate is:

```math
\hat{O}_{\mu q} = \sum_{k, n_0, n} X^{n_0 n *}_k\, f^\dagger_{k, n_0}\, f_{k+q, n} - \sum_{k, n_0, n} Y^{n_0 n *}_k\, f^\dagger_{k-q, n}\, f_{k, n_0}
```

> **Comparison with TDA**: The TDA operator is recovered by setting $Y = 0$.

### 1.2 Orthonormalization

We require orthonormality among different RPA excited states. The correct RPA normalization uses the **commutator** expectation value on the HF ground state (the quasi-boson approximation):

```math
\langle G|\, [\hat{O}_{\mu'q},\, \hat{O}^\dagger_{\mu q}]\, |G\rangle = \delta_{\mu\mu'}
```

This treats each particle-hole pair as an approximate boson: $[b, b^\dagger] \approx 1$.

Expanding $[\hat{O}_{\mu'}, \hat{O}^\dagger_\mu]$ using the bilinear commutator identity $[f^\dagger_a f_b, f^\dagger_c f_d] = \delta_{bc}\, f^\dagger_a f_d - \delta_{ad}\, f^\dagger_c f_b$:

**Forward-forward** ($X^*X$ terms): $[f^\dagger_{k',n_0'} f_{k'+q,n'},\; f^\dagger_{k+q,n} f_{k,n_0}]$ gives two bilinears. Taking $\langle G|\cdots|G\rangle$:

- $\delta_{k'+q,k+q}\,\delta_{n'n}\,\langle G| f^\dagger_{k',n_0'} f_{k,n_0} |G\rangle = \delta_{k'k}\,\delta_{n'n}\,\bar{n}_{n_0}(k)\,\delta_{n_0'n_0} = +\delta_{k'k}\,\delta_{n_0'n_0}\,\delta_{n'n}$
- $-\delta_{k'k}\,\delta_{n_0'n_0}\,\langle G| f^\dagger_{k+q,n} f_{k'+q,n'} |G\rangle = -\delta_{k'k}\,\delta_{n_0'n_0}\,\bar{n}_n(k+q)\,\delta_{nn'} = 0$ (unocc)

$$\Rightarrow\; +\sum_{k,n_0,n} X^{n_0n*}_k(\mu')\, X^{n_0n}_k(\mu)$$

**Backward-backward** ($Y^*Y$ terms): from $\hat{O}_{\mu'}$ the backward part is $-Y^*_{\mu'}\, f^\dagger_{k'-q,n'}\, f_{k',n_0'}$ and from $\hat{O}^\dagger_\mu$ the backward part is $-Y_\mu\, f^\dagger_{k,n_0}\, f_{k-q,n}$. The commutator is $[f^\dagger_{k'-q,n'} f_{k',n_0'},\; f^\dagger_{k,n_0} f_{k-q,n}]$, giving:

- $\delta_{k',k}\,\delta_{n_0',n_0}\,\langle G| f^\dagger_{k'-q,n'} f_{k-q,n} |G\rangle = \delta_{k'k}\,\delta_{n_0'n_0}\,\delta_{k'-q,k-q}\,\delta_{n'n}\,\bar{n}_{n'}(k-q) = 0$ (unocc)
- $-\delta_{k'-q,k-q}\,\delta_{n'n}\,\langle G| f^\dagger_{k,n_0} f_{k',n_0'} |G\rangle = -\delta_{k'k}\,\delta_{n'n}\,\delta_{k,k'}\,\delta_{n_0,n_0'}\,\bar{n}_{n_0}(k) = -\delta_{k'k}\,\delta_{n_0'n_0}\,\delta_{n'n}$

With the prefactor $(-Y^*_{\mu'})(-Y_\mu) = Y^*Y$:

$$\Rightarrow\; -\sum_{k,n_0,n} Y^{n_0n*}_k(\mu')\, Y^{n_0n}_k(\mu)$$

**Cross terms** ($X^*Y$ and $Y^*X$): the bilinear commutators produce terms requiring $\delta_{n_0',n}$ (occ = unocc) or $\delta_{n',n_0}$ (unocc = occ) — **vanish**.

Combining:

```math
\langle G|\, [\hat{O}_{\mu'q},\, \hat{O}^\dagger_{\mu q}]\, |G\rangle = \sum_{k, n_0, n} \left[ X^{n_0 n *}_k(\mu')\, X^{n_0 n}_k(\mu) - Y^{n_0 n *}_k(\mu')\, Y^{n_0 n}_k(\mu) \right]
```

Therefore the **RPA orthonormalization condition** is:

```math
\boxed{\sum_{k, n_0, n} \left[ X^{n_0 n *}_k(\mu')\, X^{n_0 n}_k(\mu) - Y^{n_0 n *}_k(\mu')\, Y^{n_0 n}_k(\mu) \right] = \delta_{\mu\mu'}}
```

In matrix notation: $\mathbf{X}^\dagger \mathbf{X} - \mathbf{Y}^\dagger \mathbf{Y} = \mathbf{I}$. This is a **symplectic** (indefinite) metric, in contrast to the TDA's Euclidean metric $\boldsymbol{\psi}^\dagger \boldsymbol{\psi} = \mathbf{I}$.

> **Physical meaning**: the symplectic norm reflects the bosonic commutation relation of the composite particle-hole operators. In the RPA ground state, $\langle Y^\dagger Y \rangle$ measures the depletion of the Fermi sea due to quantum fluctuations, and the norm $X^\dagger X - Y^\dagger Y = 1$ ensures unit probability after accounting for both forward and backward amplitudes.

## 2. RPA Eigenvalue Problem: Derivation from Equation of Motion

### Shorthand

We define shorthand for the bilinear operators appearing in $\hat{O}$ and $\hat{O}^\dagger$:

**From $\hat{O}^\dagger_{\mu q}$** (the excitation operator, summed over $p, n_0', n'$):
- $Q_f \equiv f^\dagger_{p+q,n'}\, f_{p,n_0'}$ — forward bilinear (coefficient $X^{n_0'n'}_p$)
- $Q_b \equiv f^\dagger_{p,n_0'}\, f_{p-q,n'}$ — backward bilinear (coefficient $-Y^{n_0'n'}_p$)

**From $\hat{O}_{\mu q}$** (the de-excitation operator, summed over $k, n_0, n$):
- $P_f \equiv f^\dagger_{k,n_0}\, f_{k+q,n}$ — forward bilinear (coefficient $X^{n_0 n *}_k$)
- $P_b \equiv f^\dagger_{k-q,n}\, f_{k,n_0}$ — backward bilinear (coefficient $-Y^{n_0 n *}_k$)

where $n_0, n_0' \in$ occ; $\quad n, n' \in$ unocc.

> **Note**: $P_f$ and $P_b$ are the bilinears appearing in $\hat{O}_{\mu q}$, and serve as the "probe" operators in the equation-of-motion method. They arise from varying $\hat{O}$ with respect to $X^*_{k,n_0,n}$ and $Y^*_{k,n_0,n}$ respectively. Note that $P_b \neq P_f^\dagger$; instead $P_f^\dagger = f^\dagger_{k+q,n}\,f_{k,n_0} = Q_f|_{p=k}$.

### 2.1 Equation of Motion

The RPA eigenvalue equation is derived from the equation-of-motion condition:

$$\langle G|\, [\delta\hat{O},\; [H,\; \hat{O}^\dagger_{\mu q}]]\, |G\rangle = \varepsilon_\mu\, \langle G|\, [\delta\hat{O},\; \hat{O}^\dagger_{\mu q}]\, |G\rangle$$

Since $\hat{O}$ is linear in $X^*$ and $Y^*$, varying with respect to each coefficient yields two independent equations:

Eq. (I) — variation w.r.t. $X^{n_0 n *}_k$ (forward probe $P_f$):

$$\langle G| \big[P_f,\; [H, \hat{O}^\dagger_{\mu q}]\big] |G\rangle = \varepsilon_\mu \langle G| \big[P_f,\; \hat{O}^\dagger_{\mu q}\big] |G\rangle$$

Eq. (II) — variation w.r.t. $Y^{n_0 n *}_k$ (backward probe $P_b$, with an extra $-$ from the $-Y^*$ coefficient):

$$-\langle G| \big[P_b,\; [H, \hat{O}^\dagger_{\mu q}]\big] |G\rangle = -\varepsilon_\mu \langle G| \big[P_b,\; \hat{O}^\dagger_{\mu q}\big] |G\rangle$$

which simplifies to:

$$\langle G| \big[P_b,\; [H, \hat{O}^\dagger_{\mu q}]\big] |G\rangle = \varepsilon_\mu \langle G| \big[P_b,\; \hat{O}^\dagger_{\mu q}\big] |G\rangle$$

---

### 2.2 Step 1: Right-Hand Sides

Substitute $\hat{O}^\dagger = \sum X \cdot Q_f - \sum Y \cdot Q_b$ and use linearity. We need the 4 commutator expectation values $\langle G|[\text{probe},\;\text{bilinear in } \hat{O}^\dagger]|G\rangle$.

#### SNEG results

$$\langle G|[P_f, Q_f]|G\rangle: \quad \delta_{k,p}\,\delta_{n,n'}\,\delta_{n_0,n_0'}\,(\theta(-k,n_0) - \theta(-k-q,n))$$

$$\langle G|[P_f, Q_b]|G\rangle: \quad \delta_{k,p-q}\,\delta_{n,n_0'}\,\delta_{n_0,n'}\,(\theta(-k,n_0) - \theta(-p,n_0'))$$

$$\langle G|[P_b, Q_f]|G\rangle: \quad \delta_{k,p+q}\,\delta_{n,n_0'}\,\delta_{n_0,n'}\,(\theta(-k+q,n) - \theta(-p-q,n'))$$

$$\langle G|[P_b, Q_b]|G\rangle: \quad \delta_{k,p}\,\delta_{n,n'}\,\delta_{n_0,n_0'}\,(\theta(-k+q,n) - \theta(-k,n_0))$$

#### Applying occupation constraints

**$\langle G|[P_f, Q_f]|G\rangle$**: $\theta(-k,n_0) - \theta(-k-q,n) = 1 - 0 = +1$

$$\Rightarrow +\delta_{kp}\,\delta_{nn'}\,\delta_{n_0 n_0'}$$

**$\langle G|[P_f, Q_b]|G\rangle$**: requires $\delta_{n,n_0'}$ (unocc = occ) **and** $\delta_{n_0,n'}$ (occ = unocc) → **vanishes**

**$\langle G|[P_b, Q_f]|G\rangle$**: requires $\delta_{n,n_0'}$ (unocc = occ) **and** $\delta_{n_0,n'}$ (occ = unocc) → **vanishes**

**$\langle G|[P_b, Q_b]|G\rangle$**: $\theta(-k+q,n) - \theta(-k,n_0) = 0 - 1 = -1$

$$\Rightarrow -\delta_{kp}\,\delta_{nn'}\,\delta_{n_0 n_0'}$$

#### Assembling the RHS

**RHS of Eq(I)**: only $[P_f, Q_f]$ survives (coefficient of $X$), $[P_f, Q_b]$ vanishes (coefficient of $-Y$):

$$\text{RHS(I)} = \varepsilon_\mu \sum_{p,n_0',n'} X^{n_0'n'}_{p} \cdot (+\delta_{kp}\delta_{nn'}\delta_{n_0 n_0'}) = \varepsilon_\mu\, X^{n_0 n}_{k}$$

**RHS of Eq(II)**: $[P_b, Q_f]$ vanishes (coefficient of $X$), $[P_b, Q_b]$ survives (coefficient of $-Y$):

$$\text{RHS(II)} = \varepsilon_\mu \sum_{p,n_0',n'} (-Y^{n_0'n'}_{p}) \cdot (-\delta_{kp}\delta_{nn'}\delta_{n_0 n_0'}) = +\varepsilon_\mu\, Y^{n_0 n}_{k}$$

**Summary**:

| Equation | RHS |
|---|---|
| Eq(I), forward probe | $\varepsilon_\mu\, X^{n_0 n}_{k}$ |
| Eq(II), backward probe | $+\varepsilon_\mu\, Y^{n_0 n}_{k}$ |

---

### 2.3 Step 2: Left-Hand Sides

Substitute $\hat{O}^\dagger = \sum X \cdot Q_f - \sum Y \cdot Q_b$ and use linearity. Define the 4 matrix elements (keeping $H$ abstract):

```math
\begin{aligned}
\mathcal{A}^{n_0 n,\,n_0'n'}_{kp} &\equiv \langle G|\,[P_f,\; [H,\; Q_f]]\,|G\rangle \\
\mathcal{B}^{n_0 n,\,n_0'n'}_{kp} &\equiv -\langle G|\,[P_f,\; [H,\; Q_b]]\,|G\rangle \\
\mathcal{C}^{n_0 n,\,n_0'n'}_{kp} &\equiv \langle G|\,[P_b,\; [H,\; Q_f]]\,|G\rangle \\
\mathcal{D}^{n_0 n,\,n_0'n'}_{kp} &\equiv -\langle G|\,[P_b,\; [H,\; Q_b]]\,|G\rangle
\end{aligned}
```

Note the minus signs in $\mathcal{B}$ and $\mathcal{D}$: they absorb the $(-Y)$ coefficient in $\hat{O}^\dagger$, so that the LHS can be written uniformly in terms of $X$ and $Y$ (without extra minus signs).

Then:

**LHS of Eq(I)**:

$$\text{LHS(I)} = \sum_{p,n_0',n'} \mathcal{A}_{kp}\, X_{p} + \sum_{p,n_0',n'} \mathcal{B}_{kp}\, Y_{p}$$

**LHS of Eq(II)**:

$$\text{LHS(II)} = \sum_{p,n_0',n'} \mathcal{C}_{kp}\, X_{p} + \sum_{p,n_0',n'} \mathcal{D}_{kp}\, Y_{p}$$

### 2.4 The RPA Eigenvalue Problem

Combining the LHS and RHS from the two equations of motion, we obtain the generalized eigenvalue problem:

```math
\boxed{
\begin{pmatrix} \mathcal{A}(q) & \mathcal{B}(q) \\ \mathcal{C}(q) & \mathcal{D}(q) \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \varepsilon_\mu
\begin{pmatrix} X \\ Y \end{pmatrix}
}
```

> **Comparison with TDA**: Setting $Y = 0$ and $\mathcal{B} = 0$ reduces to the TDA eigenvalue problem $\mathcal{A}\, X = \varepsilon_\mu\, X$. The off-diagonal blocks couple the forward and backward particle-hole channels, capturing the ground-state correlation effects absent in TDA.

> **Note on $\mathcal{A}$**: The matrix $\mathcal{A}$ is identical to the TDA effective Hamiltonian $\mathcal{H}^{\text{eff}}$ derived in the [TDA document](particle_hole.md), since the forward probe $P_f$ and forward bilinear $Q_f$ are the same as in TDA.

### 2.5 Symmetry Relation: $\mathcal{C}(q) = -\mathcal{B}(-q)^*$

Before computing each matrix explicitly, we establish the relation between $\mathcal{B}$ and $\mathcal{C}$ by expanding them with full operator indices.

#### Expanding $\mathcal{B}$ and $\mathcal{C}$

**$\mathcal{B}^{n_0 n,\,n_0'n'}_{kp}(q) = -\langle G|[P_f, [H, Q_b]]|G\rangle$:**

```math
\begin{aligned}
\mathcal{B}^{n_0 n,\,n_0'n'}_{kp}(q) =\; & -\langle G|\,f^\dagger_{k,n_0}\, f_{k+q,n}\; H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\,|G\rangle \\
& + \langle G|\,f^\dagger_{k,n_0}\, f_{k+q,n}\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; f^\dagger_{k,n_0}\, f_{k+q,n}\,|G\rangle \\
& - \langle G|\,f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\; f^\dagger_{k,n_0}\, f_{k+q,n}\,|G\rangle
\end{aligned}
```

**$\mathcal{C}^{n_0 n,\,n_0'n'}_{kp}(q) = \langle G|[P_b, [H, Q_f]]|G\rangle$:**

```math
\begin{aligned}
\mathcal{C}^{n_0 n,\,n_0'n'}_{kp}(q) =\; & +\langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\,|G\rangle \\
& - \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\,|G\rangle \\
& - \langle G|\,H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& + \langle G|\,f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle
\end{aligned}
```

#### Complex conjugate of $\mathcal{B}$

Using $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$, $H^\dagger = H$, $(ABC\cdots)^\dagger = \cdots C^\dagger B^\dagger A^\dagger$:

```math
\begin{aligned}
\mathcal{B}^{n_0 n,\,n_0'n'\,*}_{kp}(q) =\; & -\langle G|\,f^\dagger_{p-q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k+q,n}\, f_{k,n_0}\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\; f^\dagger_{k+q,n}\, f_{k,n_0}\,|G\rangle \\
& + \langle G|\,f^\dagger_{k+q,n}\, f_{k,n_0}\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\; H\,|G\rangle \\
& - \langle G|\,f^\dagger_{k+q,n}\, f_{k,n_0}\; H\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\,|G\rangle
\end{aligned}
```

#### Substituting $q \to -q$

```math
\begin{aligned}
\mathcal{B}^{n_0 n,\,n_0'n'\,*}_{kp}(-q) =\; & -\langle G|\,f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& + \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\,|G\rangle \\
& - \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\,|G\rangle
\end{aligned}
```

#### Comparing $-\mathcal{B}^*(-q)$ with $\mathcal{C}(q)$

Negating and reordering the four terms:

```math
\begin{aligned}
-\mathcal{B}^{n_0 n,\,n_0'n'\,*}_{kp}(-q) =\; & \underbrace{+\langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\,|G\rangle}_{\text{from 4th term}} \\
& \underbrace{- \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\,|G\rangle}_{\text{from 3rd term}} \\
& \underbrace{- \langle G|\,H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle}_{\text{from 2nd term}} \\
& \underbrace{+ \langle G|\,f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle}_{\text{from 1st term}}
\end{aligned}
```

This matches $\mathcal{C}^{n_0 n,\,n_0'n'}_{kp}(q)$ **term by term**.

Therefore:

```math
\boxed{\mathcal{C}(q) = -\mathcal{B}(-q)^*}
```

This relation uses only $H = H^\dagger$ and $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$. It does **not** require $|G\rangle$ to be an eigenstate of $H$, and holds for any Hermitian Hamiltonian.

### 2.6 Symmetry Relation: $\mathcal{D}(q) = -\mathcal{A}(-q)^*$

We establish the relation between $\mathcal{A}$ and $\mathcal{D}$ by the same strategy: expand with full indices, then apply complex conjugation and $q \to -q$.

#### Expanding $\mathcal{A}$ and $\mathcal{D}$

**$\mathcal{A}^{n_0 n,\,n_0'n'}_{kp}(q) = \langle G|[P_f, [H, Q_f]]|G\rangle$:**

```math
\begin{aligned}
\mathcal{A}^{n_0 n,\,n_0'n'}_{kp}(q) =\; & +\langle G|\,f^\dagger_{k,n_0}\, f_{k+q,n}\; H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\,|G\rangle \\
& - \langle G|\,f^\dagger_{k,n_0}\, f_{k+q,n}\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\,|G\rangle \\
& - \langle G|\,H\; f^\dagger_{p+q,n'}\, f_{p,n_0'}\; f^\dagger_{k,n_0}\, f_{k+q,n}\,|G\rangle \\
& + \langle G|\,f^\dagger_{p+q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k,n_0}\, f_{k+q,n}\,|G\rangle
\end{aligned}
```

**$\mathcal{D}^{n_0 n,\,n_0'n'}_{kp}(q) = -\langle G|[P_b, [H, Q_b]]|G\rangle$:**

```math
\begin{aligned}
\mathcal{D}^{n_0 n,\,n_0'n'}_{kp}(q) =\; & -\langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\,|G\rangle \\
& + \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& - \langle G|\,f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle
\end{aligned}
```

#### Substituting $q \to -q$ in $\mathcal{A}$

```math
\begin{aligned}
\mathcal{A}^{n_0 n,\,n_0'n'}_{kp}(-q) =\; & +\langle G|\,f^\dagger_{k,n_0}\, f_{k-q,n}\; H\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\,|G\rangle \\
& - \langle G|\,f^\dagger_{k,n_0}\, f_{k-q,n}\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\; H\,|G\rangle \\
& - \langle G|\,H\; f^\dagger_{p-q,n'}\, f_{p,n_0'}\; f^\dagger_{k,n_0}\, f_{k-q,n}\,|G\rangle \\
& + \langle G|\,f^\dagger_{p-q,n'}\, f_{p,n_0'}\; H\; f^\dagger_{k,n_0}\, f_{k-q,n}\,|G\rangle
\end{aligned}
```

#### Complex conjugate of $\mathcal{A}(-q)$

Using $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$, $H^\dagger = H$, $(ABC\cdots)^\dagger = \cdots C^\dagger B^\dagger A^\dagger$:

```math
\begin{aligned}
\mathcal{A}^{n_0 n,\,n_0'n'\,*}_{kp}(-q) =\; & +\langle G|\,f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& - \langle G|\,H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle \\
& - \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\,|G\rangle \\
& + \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\,|G\rangle
\end{aligned}
```

#### Comparing $-\mathcal{A}(-q)^*$ with $\mathcal{D}(q)$

Negating and reordering the four terms:

```math
\begin{aligned}
-\mathcal{A}^{n_0 n,\,n_0'n'\,*}_{kp}(-q) =\; & \underbrace{-\langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\,|G\rangle}_{\text{from 4th term}} \\
& \underbrace{+ \langle G|\,f^\dagger_{k-q,n}\, f_{k,n_0}\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\,|G\rangle}_{\text{from 3rd term}} \\
& \underbrace{+ \langle G|\,H\; f^\dagger_{p,n_0'}\, f_{p-q,n'}\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle}_{\text{from 2nd term}} \\
& \underbrace{- \langle G|\,f^\dagger_{p,n_0'}\, f_{p-q,n'}\; H\; f^\dagger_{k-q,n}\, f_{k,n_0}\,|G\rangle}_{\text{from 1st term}}
\end{aligned}
```

This matches $\mathcal{D}^{n_0 n,\,n_0'n'}_{kp}(q)$ **term by term**.

Therefore:

```math
\boxed{\mathcal{D}(q) = -\mathcal{A}(-q)^*}
```

This relation, like $\mathcal{C}(q) = -\mathcal{B}(-q)^*$, uses only $H = H^\dagger$ and $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$. It does **not** require $|G\rangle$ to be an eigenstate of $H$, and holds for any Hermitian Hamiltonian.

### 2.7 The RPA Matrix in Terms of $\mathcal{A}$ and $\mathcal{B}$

Using the two symmetry relations $\mathcal{C}(q) = -\mathcal{B}(-q)^*$ and $\mathcal{D}(q) = -\mathcal{A}(-q)^*$, the RPA eigenvalue problem §2.4 becomes:

```math
\boxed{
\begin{pmatrix} \mathcal{A}(q) & \mathcal{B}(q) \\ -\mathcal{B}(-q)^* & -\mathcal{A}(-q)^* \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \varepsilon_\mu
\begin{pmatrix} X \\ Y \end{pmatrix}
}
```

Only $\mathcal{A}(q)$ and $\mathcal{B}(q)$ need to be computed explicitly. The matrix $\mathcal{A}$ is the TDA effective Hamiltonian derived in the [TDA document](particle_hole.md); the derivation of $\mathcal{B}$ follows in §3.

---

## 3. Derivation of the $\mathcal{B}$ Matrix

The $\mathcal{B}$ matrix couples the forward and backward particle-hole channels and is defined via the double commutator evaluated in the Hartree-Fock ground state:

```math
\mathcal{B}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) \equiv -\langle G|\,[P_f,\; [H,\; Q_b]]\,|G\rangle
```

where the forward probe is $P_f = f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}+\mathbf{q},n}$ and the backward bilinear is $Q_b = f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}-\mathbf{q},n'}$. Here $n_0, n_0' \in$ occ and $n, n' \in$ unocc.

### 3.1 The One-Body Parts Vanish

We decompose the Hamiltonian as $H = H_0 + H_{\text{SE}} + :H_{\text{int}}:$ (cf. [TDA document §3.1](particle_hole.md)). Since $H_0$ and $H_{\text{SE}}$ are one-body operators, $[H_{\text{1-body}}, Q_b]$ produces bilinears, and the expectation value $\langle G|[P_f, \text{bilinear}]|G\rangle$ reduces to products of Kronecker deltas. The non-vanishing contractions would require occupied indices to equal unoccupied indices (e.g., $n_0 = n'$ or $n = n_0'$), which is strictly impossible. Thus, **all one-body contributions vanish identically**:

```math
\mathcal{B}_0 = \mathcal{B}_{\text{SE}} = 0
```

> **Contrast with $\mathcal{A}$**: In the TDA matrix $\mathcal{A}$, the one-body parts contribute the diagonal mean-field energy $\delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k})$. This is possible because $P_f$ and $Q_f$ share the same occ/unocc structure. In $\mathcal{B}$, $P_f$ and $Q_b$ have reversed occ/unocc roles, making all one-body contractions vanish.

### 3.2 Normal-Ordered Two-Body Part: Expanding the Inner Commutator

The entire contribution comes from the normal-ordered interaction:

```math
\mathcal{B} = -\langle G|\,[P_f,\; [:H_{\text{int}}:,\; Q_b]]\,|G\rangle
```

Using the SNEG-verified commutator identity for $[f^\dagger_1 f^\dagger_3 f_2 f_4,\; f^\dagger_5 f_6]$ with $5 = (\mathbf{p}, n_0')$ and $6 = (\mathbf{p}-\mathbf{q}, n')$, we obtain four terms:

| Term | Commutator delta | Sign | Remaining normal-ordered string |
|---|---|---|---|
| (A) | $\delta_{\mathbf{p}_1,\mathbf{p}-\mathbf{q}}\,\delta_{m_1,n'}$ | $+$ | $f^\dagger_{\mathbf{p}_3,m_3}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$ |
| (B) | $\delta_{\mathbf{p}_2,\mathbf{p}}\,\delta_{m_2,n_0'}$ | $-$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}\, f_{\mathbf{p}-\mathbf{q},n'}$ |
| (C) | $\delta_{\mathbf{p}_3,\mathbf{p}-\mathbf{q}}\,\delta_{m_3,n'}$ | $-$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$ |
| (D) | $\delta_{\mathbf{p}_4,\mathbf{p}}\,\delta_{m_4,n_0'}$ | $+$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}-\mathbf{q},n'}$ |

### 3.3 Vanishing of Terms (B) and (D)

Terms (B) and (D) contain the annihilation operator $f_{\mathbf{p}-\mathbf{q},n'}$ at the far right of the normal-ordered string. Since $n' \in \text{unocc}$, its occupation number vanishes: $\bar{n}_{n'}(\mathbf{p}-\mathbf{q}) = 0$. In the subsequent Wick contraction with $\langle G|\, P_f \cdots |G\rangle$, $f_{\mathbf{p}-\mathbf{q},n'}$ can only pair with a creation operator to its left, producing a contraction of the form $\langle G| f^\dagger_X\, f_{\mathbf{p}-\mathbf{q},n'} |G\rangle = \delta_X \cdot \bar{n}_{n'}(\mathbf{p}-\mathbf{q}) = 0$. Therefore, **Terms (B) and (D) vanish identically**.

> **Note**: This is the key simplification that reduces the $\mathcal{B}$ computation from four terms to two.

### 3.4 Wick Contractions for Terms (A) and (C)

For the surviving Terms (A) and (C), we sandwich with $-\langle G|\, f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}+\mathbf{q},n} \cdots |G\rangle$ and evaluate the six-operator expectation values via Wick's theorem.

The external unoccupied annihilator $f_{\mathbf{k}+\mathbf{q},n}$ ($n \in$ unocc) must pair with one of the two internal $f^\dagger$'s, giving two sub-contractions per term. The overall prefactor is $(-1)_{\mathcal{B}} \times (-1/N)_{:H_{\text{int}}:} = +1/N$.

#### Term (A): $+\delta_{\mathbf{p}_1,\mathbf{p}-\mathbf{q}}\,\delta_{m_1,n'}$

The commutator delta fixes $\mathbf{p}_1 = \mathbf{p}-\mathbf{q}$, $m_1 = n'$. The remaining string is $f^\dagger_{\mathbf{p}_3,m_3}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$.

The external $f_{\mathbf{k}+\mathbf{q},n}$ must pair with $f^\dagger_{\mathbf{p}_3,m_3}$, fixing $\mathbf{p}_3 = \mathbf{k}+\mathbf{q}$, $m_3 = n$. Two sub-contractions remain for the external $f^\dagger_{\mathbf{k},n_0}$ and $f^\dagger_{\mathbf{p},n_0'}$ pairing with $f_{\mathbf{p}_2,m_2}$ and $f_{\mathbf{p}_4,m_4}$:

**A-1 (direct)**: Pair $f^\dagger_{\mathbf{k},n_0}$ with $f_{\mathbf{p}_2,m_2}$ and $f^\dagger_{\mathbf{p},n_0'}$ with $f_{\mathbf{p}_4,m_4}$.

Resolving: $\mathbf{p}_2 = \mathbf{k}$, $m_2 = n_0$; $\mathbf{p}_4 = \mathbf{p}$, $m_4 = n_0'$. Momentum conservation: $\mathbf{p}_4 = \mathbf{p}_1 + \mathbf{p}_3 - \mathbf{p}_2 = (\mathbf{p}-\mathbf{q}) + (\mathbf{k}+\mathbf{q}) - \mathbf{k} = \mathbf{p}$ ✓

Wick sign: $-1$ (one crossing). Total sign: $(+1/N) \times (+1)_{\text{comm}} \times (-1)_{\text{Wick}} = -1/N$.

$$\Rightarrow -\frac{1}{N}\,\widetilde{V}_{n',\, n_0,\, n,\, n_0'}(\mathbf{p}-\mathbf{q},\, \mathbf{k},\, \mathbf{k}+\mathbf{q})$$

**A-2 (exchange)**: Pair $f^\dagger_{\mathbf{k},n_0}$ with $f_{\mathbf{p}_4,m_4}$ and $f^\dagger_{\mathbf{p},n_0'}$ with $f_{\mathbf{p}_2,m_2}$.

Resolving: $\mathbf{p}_4 = \mathbf{k}$, $m_4 = n_0$; $\mathbf{p}_2 = \mathbf{p}$, $m_2 = n_0'$. Momentum conservation: $\mathbf{p}_4 = (\mathbf{p}-\mathbf{q}) + (\mathbf{k}+\mathbf{q}) - \mathbf{p} = \mathbf{k}$ ✓

Wick sign: $+1$ (no crossing). Total sign: $(+1/N) \times (+1)_{\text{comm}} \times (+1)_{\text{Wick}} = +1/N$.

$$\Rightarrow +\frac{1}{N}\,\widetilde{V}_{n',\, n_0',\, n,\, n_0}(\mathbf{p}-\mathbf{q},\, \mathbf{p},\, \mathbf{k}+\mathbf{q})$$

#### Term (C): $-\delta_{\mathbf{p}_3,\mathbf{p}-\mathbf{q}}\,\delta_{m_3,n'}$

The commutator delta fixes $\mathbf{p}_3 = \mathbf{p}-\mathbf{q}$, $m_3 = n'$, with an overall sign $-1$. The remaining string is $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$.

The external $f_{\mathbf{k}+\mathbf{q},n}$ must pair with $f^\dagger_{\mathbf{p}_1,m_1}$, fixing $\mathbf{p}_1 = \mathbf{k}+\mathbf{q}$, $m_1 = n$.

**C-1 (exchange)**: Pair $f^\dagger_{\mathbf{k},n_0}$ with $f_{\mathbf{p}_2,m_2}$ and $f^\dagger_{\mathbf{p},n_0'}$ with $f_{\mathbf{p}_4,m_4}$.

Resolving: $\mathbf{p}_2 = \mathbf{k}$, $m_2 = n_0$; $\mathbf{p}_4 = \mathbf{p}$, $m_4 = n_0'$. Momentum conservation: $\mathbf{p}_4 = (\mathbf{k}+\mathbf{q}) + (\mathbf{p}-\mathbf{q}) - \mathbf{k} = \mathbf{p}$ ✓

Wick sign: $-1$. Total sign: $(+1/N) \times (-1)_{\text{comm}} \times (-1)_{\text{Wick}} = +1/N$.

$$\Rightarrow +\frac{1}{N}\,\widetilde{V}_{n,\, n_0,\, n',\, n_0'}(\mathbf{k}+\mathbf{q},\, \mathbf{k},\, \mathbf{p}-\mathbf{q})$$

**C-2 (direct)**: Pair $f^\dagger_{\mathbf{k},n_0}$ with $f_{\mathbf{p}_4,m_4}$ and $f^\dagger_{\mathbf{p},n_0'}$ with $f_{\mathbf{p}_2,m_2}$.

Resolving: $\mathbf{p}_4 = \mathbf{k}$, $m_4 = n_0$; $\mathbf{p}_2 = \mathbf{p}$, $m_2 = n_0'$. Momentum conservation: $\mathbf{p}_4 = (\mathbf{k}+\mathbf{q}) + (\mathbf{p}-\mathbf{q}) - \mathbf{p} = \mathbf{k}$ ✓

Wick sign: $+1$. Total sign: $(+1/N) \times (-1)_{\text{comm}} \times (+1)_{\text{Wick}} = -1/N$.

$$\Rightarrow -\frac{1}{N}\,\widetilde{V}_{n,\, n_0',\, n',\, n_0}(\mathbf{k}+\mathbf{q},\, \mathbf{p},\, \mathbf{p}-\mathbf{q})$$

#### Topology Classification

| Sub-term | $\widetilde{V}$ band indices | $\widetilde{V}$ momenta $(\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_3)$ | Sign | Topology |
|---|---|---|---|---|
| A-1 | $(n', n_0, n, n_0')$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$ | $-1/N$ | Direct |
| C-2 | $(n, n_0', n', n_0)$ | $(\mathbf{k}+\mathbf{q},\; \mathbf{p},\; \mathbf{p}-\mathbf{q})$ | $-1/N$ | Direct |
| A-2 | $(n', n_0', n, n_0)$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{p},\; \mathbf{k}+\mathbf{q})$ | $+1/N$ | Exchange |
| C-1 | $(n, n_0, n', n_0')$ | $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p}-\mathbf{q})$ | $+1/N$ | Exchange |

Each channel receives two contributions with different $\widetilde{V}$ momentum arguments. The signs match the $\mathcal{A}$ kernel: **direct $-1/N$, exchange $+1/N$**.

### 3.5 Explicit Kernel Expressions in Orbital Basis

Transforming back to the orbital basis using $\widetilde{V}_{m_1 m_2 m_3 m_4}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3) = \sum_{abcd} U^*_{am_1}(\mathbf{p}_1)\, U_{bm_2}(\mathbf{p}_2)\, U^*_{cm_3}(\mathbf{p}_3)\, U_{dm_4}(\mathbf{p}_4)\; \widetilde{V}^{abcd}(\mathbf{p}_1,\mathbf{p}_2,\mathbf{p}_3)$:

**(A) Direct channel** $\mathcal{B}^{\text{d}}$:

```math
\boxed{
\begin{aligned}
\mathcal{B}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = -\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an'}(\mathbf{p}-\mathbf{q})\, U_{bn_0}(\mathbf{k})\, U^*_{cn}(\mathbf{k}+\mathbf{q})\, U_{dn_0'}(\mathbf{p})\; \widetilde{V}^{abcd}(\mathbf{p}-\mathbf{q},\, \mathbf{k},\, \mathbf{k}+\mathbf{q}) \\
+\; & U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{bn_0'}(\mathbf{p})\, U^*_{cn'}(\mathbf{p}-\mathbf{q})\, U_{dn_0}(\mathbf{k})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q},\, \mathbf{p},\, \mathbf{p}-\mathbf{q})
\Big]
\end{aligned}
}
```

**(B) Exchange channel** $\mathcal{B}^{\text{x}}$:

```math
\boxed{
\begin{aligned}
\mathcal{B}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = +\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an'}(\mathbf{p}-\mathbf{q})\, U_{bn_0'}(\mathbf{p})\, U^*_{cn}(\mathbf{k}+\mathbf{q})\, U_{dn_0}(\mathbf{k})\; \widetilde{V}^{abcd}(\mathbf{p}-\mathbf{q},\, \mathbf{p},\, \mathbf{k}+\mathbf{q}) \\
+\; & U^*_{an}(\mathbf{k}+\mathbf{q})\, U_{bn_0}(\mathbf{k})\, U^*_{cn'}(\mathbf{p}-\mathbf{q})\, U_{dn_0'}(\mathbf{p})\; \widetilde{V}^{abcd}(\mathbf{k}+\mathbf{q},\, \mathbf{k},\, \mathbf{p}-\mathbf{q})
\Big]
\end{aligned}
}
```

> **Structural comparison with $\mathcal{A}$**: The $\mathcal{B}$ kernel has the same sign structure as $\mathcal{A}$ (direct $-1/N$, exchange $+1/N$) and the same four external momenta $\{\mathbf{k}, \mathbf{k}+\mathbf{q}, \mathbf{p}, \cdot\}$. The only difference is that wherever $\mathbf{p}+\mathbf{q}$ appears in $\mathcal{A}$, the corresponding position in $\mathcal{B}$ is replaced by $\mathbf{p}-\mathbf{q}$. This reflects the reversed momentum structure of the backward bilinear $Q_b$.

---

## 4. Hermiticity and Inter-Block Relations for Implementation

In §2.5–2.6 we derived the symmetry relations $\mathcal{C}(q) = -\mathcal{B}(-q)^*$ and $\mathcal{D}(q) = -\mathcal{A}(-q)^*$. These involve the matrices evaluated at $-q$, which may have a **different** particle-hole pair space than at $+q$ (different occupation at $\mathbf{k}+\mathbf{q}$ vs. $\mathbf{k}-\mathbf{q}$). This complicates implementation. In this section we derive **alternative relations** that express $\mathcal{C}$ purely in terms of $\mathcal{B}(q)$, and establish hermiticity properties of each block.

### 4.1 Relation $\mathcal{C}(q) = -\mathcal{B}(q)^\dagger$

**Claim:**

```math
\boxed{\mathcal{C}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p}}(q) = -\mathcal{B}^{n_0'n',\,n_0 n\,*}_{\mathbf{p}\mathbf{k}}(q)}
```

i.e., $\mathcal{C}(q) = -\mathcal{B}(q)^\dagger$.

**Proof.** We follow the same strategy as §2.5: expand with full operator indices, take the complex conjugate, and compare term by term.

#### Step 1: Expanding $\mathcal{B}^{n_0'n',\,n_0 n}_{\mathbf{p}\mathbf{k}}(q)$

The matrix element $\mathcal{B}^{n_0'n',\,n_0 n}_{\mathbf{p}\mathbf{k}}(q)$ is obtained from the definition $\mathcal{B} = -\langle G|[P_f, [H, Q_b]]|G\rangle$ by substituting $(\mathbf{k}, n_0, n) \to (\mathbf{p}, n_0', n')$ and $(\mathbf{p}, n_0', n') \to (\mathbf{k}, n_0, n)$, giving $P_f' = f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}+\mathbf{q},n'}$ and $Q_b' = f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}-\mathbf{q},n}$:

```math
\begin{aligned}
\mathcal{B}^{n_0'n',\,n_0 n}_{\mathbf{p}\mathbf{k}}(q) =\; & -\langle G|\,f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}+\mathbf{q},n'}\; H\; f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}-\mathbf{q},n}\,|G\rangle \\
& + \langle G|\,f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}+\mathbf{q},n'}\; f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}-\mathbf{q},n}\; H\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}-\mathbf{q},n}\; f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}+\mathbf{q},n'}\,|G\rangle \\
& - \langle G|\,f^\dagger_{\mathbf{k},n_0}\, f_{\mathbf{k}-\mathbf{q},n}\; H\; f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}+\mathbf{q},n'}\,|G\rangle
\end{aligned}
```

#### Step 2: Complex conjugate

Using $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$, $H^\dagger = H$, $(ABC\cdots)^\dagger = \cdots C^\dagger B^\dagger A^\dagger$:

```math
\begin{aligned}
\mathcal{B}^{n_0'n',\,n_0 n\,*}_{\mathbf{p}\mathbf{k}}(q) =\; & -\langle G|\,f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; H\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\,|G\rangle \\
& + \langle G|\,H\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\,|G\rangle \\
& + \langle G|\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; H\,|G\rangle \\
& - \langle G|\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; H\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\,|G\rangle
\end{aligned}
```

#### Step 3: Comparing $-\mathcal{B}^{n_0'n',\,n_0 n\,*}_{\mathbf{p}\mathbf{k}}(q)$ with $\mathcal{C}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p}}(q)$

Negating and reordering the four terms:

```math
\begin{aligned}
-\mathcal{B}^{n_0'n',\,n_0 n\,*}_{\mathbf{p}\mathbf{k}}(q) =\; & \underbrace{+\langle G|\,f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; H\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\,|G\rangle}_{\text{from 1st term}} \\
& \underbrace{- \langle G|\,f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; H\,|G\rangle}_{\text{from 3rd term (negated)}} \\
& \underbrace{- \langle G|\,H\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\,|G\rangle}_{\text{from 3rd term}} \\
& \underbrace{+ \langle G|\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; H\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\,|G\rangle}_{\text{from 4th term}}
\end{aligned}
```

Now recall $\mathcal{C}$ from §2.5:

```math
\begin{aligned}
\mathcal{C}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p}}(q) =\; & +\langle G|\,f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; H\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\,|G\rangle \\
& - \langle G|\,f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; H\,|G\rangle \\
& - \langle G|\,H\; f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\,|G\rangle \\
& + \langle G|\,f^\dagger_{\mathbf{p}+\mathbf{q},n'}\, f_{\mathbf{p},n_0'}\; H\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}\,|G\rangle
\end{aligned}
```

The four terms match **term by term**.

Therefore:

```math
\boxed{\mathcal{C}(q) = -\mathcal{B}(q)^\dagger}
```

This relation uses only $H = H^\dagger$ and $\langle G|X|G\rangle^* = \langle G|X^\dagger|G\rangle$. It does **not** require $|G\rangle$ to be an eigenstate of $H$, and holds for any Hermitian Hamiltonian. Unlike the §2.5 relation $\mathcal{C}(q) = -\mathcal{B}(-q)^*$, this form involves only $\mathcal{B}$ at the **same** momentum transfer $q$, making it directly usable for implementation.

### 4.2 Relation Between $\mathcal{D}(q)$ and $\mathcal{A}(-q)$

From §2.6, we have $\mathcal{D}(q) = -\mathcal{A}(-q)^*$. This involves $\mathcal{A}$ at $-q$, which uses eigensystems at $\mathbf{k}-\mathbf{q}$ and $\mathbf{p}-\mathbf{q}$ — genuinely different momentum points from the $+q$ case. No relabelling of summation momenta can avoid this: the backward bilinears $P_b$ and $Q_b$ intrinsically involve $\mathbf{k}-\mathbf{q}$ and $\mathbf{p}-\mathbf{q}$, distinguishing $\mathcal{D}$ from $\mathcal{A}$.

When the particle-hole pair space is the same at $+q$ and $-q$ (i.e., $\text{unocc}(\mathbf{k}+\mathbf{q}) = \text{unocc}(\mathbf{k}-\mathbf{q})$ for all $\mathbf{k}$), one can simply compute $\mathcal{A}(-q)$ and use $\mathcal{D}(q) = -\mathcal{A}(-q)^*$. When the pair spaces differ, $\mathcal{D}(q)$ must be computed directly from its own kernel expressions, which we derive in §5.

---

## 5. Derivation of the $\mathcal{D}$ Matrix

The $\mathcal{D}$ matrix couples the backward probe and backward bilinear:

```math
\mathcal{D}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) \equiv -\langle G|\,[P_b,\; [H,\; Q_b]]\,|G\rangle
```

where $P_b = f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}$ and $Q_b = f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}-\mathbf{q},n'}$, with $n_0, n_0' \in$ occ and $n, n' \in$ unocc.

> **Note on the overall sign**: The $\mathcal{D}$ definition carries a minus sign (§2.3), so $\mathcal{D} = -\langle G|[P_b, [H, Q_b]]|G\rangle$. Combined with the $-1/N$ from $:H_{\text{int}}:$, the overall prefactor for the two-body part is $(-1) \times (-1/N) = +1/N$ — the **same** as for $\mathcal{B}$ and **opposite** to $\mathcal{A}$ (which has $+1 \times (-1/N) = -1/N$).

### 5.1 The One-Body Part

The one-body Hamiltonian $H_0 + H_{\text{SE}}$ gives the diagonal mean-field energy. Since $P_b$ and $Q_b$ share the same occ/unocc pairing structure as $P_f$ and $Q_f$ (each pairs one occupied and one unoccupied index), the one-body contributions **do not vanish** (unlike in $\mathcal{B}$):

```math
\mathcal{D}^{n_0 n,\,n_0'n'}_{\mathbf{k}\mathbf{p},\,\text{MF}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^{n_0}_\mathbf{k} - E^n_{\mathbf{k}-\mathbf{q}}\right)
```

> **Comparison with $\mathcal{A}$**: $\mathcal{A}_{\text{MF}} = \delta(E^n_{\mathbf{k}+\mathbf{q}} - E^{n_0}_\mathbf{k})$, which is positive for particle-hole excitations. In $\mathcal{D}$, the sign is reversed and $\mathbf{q} \to -\mathbf{q}$: $\mathcal{D}_{\text{MF}} = \delta(E^{n_0}_\mathbf{k} - E^n_{\mathbf{k}-\mathbf{q}})$, which is generically **negative**.

### 5.2 Normal-Ordered Two-Body Part: Inner Commutator

The entire interaction contribution comes from $:H_{\text{int}}:$:

```math
\mathcal{D}_{\text{int}} = -\langle G|\,[P_b,\; [:H_{\text{int}}:,\; Q_b]]\,|G\rangle
```

The inner commutator $[:H_{\text{int}}:, Q_b]$ is **identical** to the one computed in §3.2 for $\mathcal{B}$, since the same $Q_b$ appears in both. The four terms from the table in §3.2 apply directly:

| Term | Commutator delta | Sign | Remaining normal-ordered string |
|---|---|---|---|
| (A) | $\delta_{\mathbf{p}_1,\mathbf{p}-\mathbf{q}}\,\delta_{m_1,n'}$ | $+$ | $f^\dagger_{\mathbf{p}_3,m_3}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$ |
| (B) | $\delta_{\mathbf{p}_2,\mathbf{p}}\,\delta_{m_2,n_0'}$ | $-$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_4,m_4}\, f_{\mathbf{p}-\mathbf{q},n'}$ |
| (C) | $\delta_{\mathbf{p}_3,\mathbf{p}-\mathbf{q}}\,\delta_{m_3,n'}$ | $-$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p},n_0'}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}_4,m_4}$ |
| (D) | $\delta_{\mathbf{p}_4,\mathbf{p}}\,\delta_{m_4,n_0'}$ | $+$ | $f^\dagger_{\mathbf{p}_1,m_1}\, f^\dagger_{\mathbf{p}_3,m_3}\, f_{\mathbf{p}_2,m_2}\, f_{\mathbf{p}-\mathbf{q},n'}$ |

### 5.3 Vanishing of Terms (B) and (D)

Terms (B) and (D) contain $f_{\mathbf{p}-\mathbf{q},n'}$ ($n' \in$ unocc) at the far right, exactly as in §3.3. By the same argument: $\bar{n}_{n'}(\mathbf{p}-\mathbf{q}) = 0$, so all Wick contractions involving this operator vanish. **Terms (B) and (D) vanish identically.**

### 5.4 Wick Contractions for Terms (A) and (C)

To evaluate $\mathcal{D}_{\text{int}} = -\langle G|\,[P_b,\; X]\,|G\rangle$ (where $X$ represents the surviving terms of the inner commutator), we expand the outer commutator:

$$\mathcal{D}_{\text{int}} = -\langle G|\, P_b X \,|G\rangle + \langle G|\, X P_b \,|G\rangle$$

Notice that $P_b = f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0}$. Because $n$ is an unoccupied band, the creation operator $f^\dagger_{\mathbf{k}-\mathbf{q},n}$ acting to the left on the bra annihilates the Hartree-Fock ground state (i.e., $\langle G| f^\dagger_{\text{unocc}} = 0$). Therefore, the first term strictly vanishes: $\langle G|\, P_b X \,|G\rangle = 0$. 

The entire $\mathcal{D}$ matrix contribution thus comes elegantly from the second term, where the minus sign naturally disappears:

$$\mathcal{D}_{\text{int}} = +\langle G|\, X\; f^\dagger_{\mathbf{k}-\mathbf{q},n}\, f_{\mathbf{k},n_0} \,|G\rangle$$

Here, the outer probe $P_b$ acts on the **right**. The overall prefactor is $(+1) \times (-1/N)_{:H_{\text{int}}:} = -1/N$. 
For a non-zero Wick contraction, the external occ annihilator $f_{\mathbf{k},n_0}$ (at the far right) must pair with an internal $f^\dagger$ (creator), and the external unocc creator $f^\dagger_{\mathbf{k}-\mathbf{q},n}$ must pair with an internal $f$ (annihilator).

#### SNEG verification

The Wick contractions were computed using the SNEG code:

```mathematica
<< sneg`sneg`;
snegfermionoperators[f];
ordering[f] = SEA;
snegfreeindexes[k, p, q, n0, n, n0p, np, p1, p2, p3, m1, m2, m3, m4];

Pb = nc[f[CR, k - q, n], f[AN, k, n0]];

(* Term A: delta(p1,p-q) delta(m1,n') applied *)
termA = nc[f[CR, p3, m3], f[CR, p, n0p], f[AN, p2, m2],
           f[AN, p - q + p3 - p2, m4]];

(* Term C: delta(p3,p-q) delta(m3,n') applied *)
termC = nc[f[CR, p1, m1], f[CR, p, n0p], f[AN, p2, m2],
           f[AN, p1 + p - q - p2, m4]];

commA = Expand[komutator[Pb, termA]];
RA = SimplifyKD[vevwick[commA]];
Print["R_A = ", RA];

commC = Expand[komutator[Pb, termC]];
RC = SimplifyKD[vevwick[commC]];
Print["R_C = ", RC];
```

Each result (R_A, R_C) contains two groups:

- **Group 1** ($\delta_{\mathbf{k},\mathbf{p}}\,\delta_{n_0,n_0'}$): self-energy terms
- **Group 2** ($\delta_{\mathbf{k},\mathbf{p}_3}\,\delta_{m_3,n_0}$ for R_A, or $\delta_{\mathbf{k},\mathbf{p}_1}\,\delta_{m_1,n_0}$ for R_C): scattering terms

#### Understanding the $\theta$ factors

Unlike $\mathcal{A}$ and $\mathcal{B}$, the scattering terms in $\mathcal{D}$ carry a **residual occupation factor**. For R_A's scattering group, SNEG gives:

```math
\theta(-\mathbf{p}, n_0')\,\big(\theta(-\mathbf{p}_3, n_0) - \theta(-\mathbf{k}+\mathbf{q}, n)\big)
```

After resolving $\mathbf{p}_3 = \mathbf{k}$, $m_3 = n_0$:
- $\theta(-\mathbf{p}, n_0') = \bar{n}_{n_0'}(\mathbf{p}) = 1$ (since $n_0' \in$ occ)
- $\theta(-\mathbf{k}, n_0) = \bar{n}_{n_0}(\mathbf{k}) = 1$ (since $n_0 \in$ occ)
- $\theta(-(\mathbf{k}-\mathbf{q}), n) = \bar{n}_n(\mathbf{k}-\mathbf{q})$ — this **cannot** be simplified in general

The factor becomes $1 - \bar{n}_n(\mathbf{k}-\mathbf{q})$, which equals:
- **1** if $n \in \text{unocc}(\mathbf{k}-\mathbf{q})$
- **0** if $n \in \text{occ}(\mathbf{k}-\mathbf{q})$

> **Contrast with $\mathcal{A}$ and $\mathcal{B}$**: In those matrices, every $\theta$ factor involving an external index is either trivially 0 or 1 (e.g., $\bar{n}_n(\mathbf{k}+\mathbf{q}) = 0$ because $n \in \text{unocc}(\mathbf{k}+\mathbf{q})$ by construction). In $\mathcal{D}$, the residual factor $1 - \bar{n}_n(\mathbf{k}-\mathbf{q})$ reflects the fact that $n$ is defined as unoccupied at $\mathbf{k}+\mathbf{q}$, but its occupation at $\mathbf{k}-\mathbf{q}$ is not constrained.

**In what follows, we assume that $n \in \text{unocc}(\mathbf{k}-\mathbf{q})$** (and similarly $n' \in \text{unocc}(\mathbf{p}-\mathbf{q})$), so that all $\theta$ factors reduce to 1. This is automatically satisfied when the particle-hole triples are defined by requiring $n$ to be unoccupied at both $\mathbf{k}+\mathbf{q}$ and $\mathbf{k}-\mathbf{q}$. When this condition does not hold (i.e., band $n$ is occupied at $\mathbf{k}-\mathbf{q}$ but unoccupied at $\mathbf{k}+\mathbf{q}$), the corresponding $\mathcal{D}$ matrix elements vanish.

#### Scattering contributions from R_A

The SNEG scattering group of R_A (after setting all $\theta = 1$) gives two sub-terms, which combined with the overall prefactor $+1/N$ yield:

**A-direct**: $\widetilde{V}_{n',\, n,\, n_0,\, n_0'}(\mathbf{p}-\mathbf{q},\, \mathbf{k}-\mathbf{q},\, \mathbf{k})$, sign $+1/N$

**A-exchange**: $\widetilde{V}_{n',\, n_0',\, n_0,\, n}(\mathbf{p}-\mathbf{q},\, \mathbf{p},\, \mathbf{k})$, sign $-1/N$

#### Scattering contributions from R_C

By the same analysis (R_C has identical structure with $\mathbf{p}_3 \to \mathbf{p}_1$), combined with the commutator sign $-1$ for Term C:

**C-direct**: $\widetilde{V}_{n_0,\, n_0',\, n',\, n}(\mathbf{k},\, \mathbf{p},\, \mathbf{p}-\mathbf{q})$, sign $+1/N$

**C-exchange**: $\widetilde{V}_{n_0,\, n,\, n',\, n_0'}(\mathbf{k},\, \mathbf{k}-\mathbf{q},\, \mathbf{p}-\mathbf{q})$, sign $-1/N$

#### Topology Classification

| Sub-term | $\widetilde{V}$ band indices | $\widetilde{V}$ momenta $(\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_3)$ | Sign | Topology |
|---|---|---|---|---|
| A-direct | $(n', n, n_0, n_0')$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{k}-\mathbf{q},\; \mathbf{k})$ | $+1/N$ | Direct |
| C-direct | $(n_0, n_0', n', n)$ | $(\mathbf{k},\; \mathbf{p},\; \mathbf{p}-\mathbf{q})$ | $+1/N$ | Direct |
| A-exchange | $(n', n_0', n_0, n)$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{p},\; \mathbf{k})$ | $-1/N$ | Exchange |
| C-exchange | $(n_0, n, n', n_0')$ | $(\mathbf{k},\; \mathbf{k}-\mathbf{q},\; \mathbf{p}-\mathbf{q})$ | $-1/N$ | Exchange |

Each channel has two contributions. The signs are **opposite** to $\mathcal{A}$ and $\mathcal{B}$: **direct $+1/N$, exchange $-1/N$**. This is consistent with the symmetry relation $\mathcal{D}(q) = -\mathcal{A}(-q)^*$: the negation inverts the sign structure.

### 5.5 Explicit Kernel Expressions in Orbital Basis

Transforming to the orbital basis:

**(A) Direct channel** $\mathcal{D}^{\text{d}}$:

```math
\boxed{
\begin{aligned}
\mathcal{D}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = +\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an'}(\mathbf{p}-\mathbf{q})\, U_{bn}(\mathbf{k}-\mathbf{q})\, U^*_{cn_0}(\mathbf{k})\, U_{dn_0'}(\mathbf{p})\; \widetilde{V}^{abcd}(\mathbf{p}-\mathbf{q},\, \mathbf{k}-\mathbf{q},\, \mathbf{k}) \\
+\; & U^*_{an_0}(\mathbf{k})\, U_{bn_0'}(\mathbf{p})\, U^*_{cn'}(\mathbf{p}-\mathbf{q})\, U_{dn}(\mathbf{k}-\mathbf{q})\; \widetilde{V}^{abcd}(\mathbf{k},\, \mathbf{p},\, \mathbf{p}-\mathbf{q})
\Big]
\end{aligned}
}
```

**(B) Exchange channel** $\mathcal{D}^{\text{x}}$:

```math
\boxed{
\begin{aligned}
\mathcal{D}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = -\frac{1}{N}\sum_{abcd} \Big[
& U^*_{an'}(\mathbf{p}-\mathbf{q})\, U_{bn_0'}(\mathbf{p})\, U^*_{cn_0}(\mathbf{k})\, U_{dn}(\mathbf{k}-\mathbf{q})\; \widetilde{V}^{abcd}(\mathbf{p}-\mathbf{q},\, \mathbf{p},\, \mathbf{k}) \\
+\; & U^*_{an_0}(\mathbf{k})\, U_{bn}(\mathbf{k}-\mathbf{q})\, U^*_{cn'}(\mathbf{p}-\mathbf{q})\, U_{dn_0'}(\mathbf{p})\; \widetilde{V}^{abcd}(\mathbf{k},\, \mathbf{k}-\mathbf{q},\, \mathbf{p}-\mathbf{q})
\Big]
\end{aligned}
}
```

> **Structural comparison with $\mathcal{A}$**: The $\mathcal{D}$ kernel has **opposite** sign structure to $\mathcal{A}$: direct $+1/N$ (vs. $-1/N$ in $\mathcal{A}$) and exchange $-1/N$ (vs. $+1/N$ in $\mathcal{A}$). This sign reversal is the manifestation of $\mathcal{D}(q) = -\mathcal{A}(-q)^*$. Comparing the momentum arguments:
>
> | | $\mathcal{A}$ momenta | $\mathcal{D}$ momenta |
> |---|---|---|
> | Direct-1 | $(\mathbf{k}+\mathbf{q},\; \mathbf{p}+\mathbf{q},\; \mathbf{p})$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{k}-\mathbf{q},\; \mathbf{k})$ |
> | Direct-2 | $(\mathbf{p},\; \mathbf{k},\; \mathbf{k}+\mathbf{q})$ | $(\mathbf{k},\; \mathbf{p},\; \mathbf{p}-\mathbf{q})$ |
> | Exchange-1 | $(\mathbf{k}+\mathbf{q},\; \mathbf{k},\; \mathbf{p})$ | $(\mathbf{p}-\mathbf{q},\; \mathbf{p},\; \mathbf{k})$ |
> | Exchange-2 | $(\mathbf{p},\; \mathbf{p}+\mathbf{q},\; \mathbf{k}+\mathbf{q})$ | $(\mathbf{k},\; \mathbf{k}-\mathbf{q},\; \mathbf{p}-\mathbf{q})$ |
>
> In $\mathcal{D}$, all $+\mathbf{q}$ shifts are replaced by $-\mathbf{q}$, and the roles of $(\mathbf{k}, n_0)$ and $(\mathbf{p}, n_0')$ are swapped relative to $\mathcal{A}$. The band indices also differ accordingly.

### 5.6 Final Result

```math
\boxed{
\mathcal{D}^{n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) = \delta_{\mathbf{k}\mathbf{p}}\,\delta_{n_0 n_0'}\,\delta_{nn'}\left(E^{n_0}_\mathbf{k} - E^n_{\mathbf{k}-\mathbf{q}}\right) + \mathcal{D}^{\text{d},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q}) + \mathcal{D}^{\text{x},\,n_0 n,\, n_0' n'}_{\mathbf{k}\mathbf{p}}(\mathbf{q})
}
```

### 5.7 Summary: Practical RPA Matrix Assembly

Combining all results from §3–5, the RPA eigenvalue problem can be assembled as:

```math
\boxed{
\begin{pmatrix} \mathcal{A}(\mathbf{q}) & \mathcal{B}(\mathbf{q}) \\ -\mathcal{B}(\mathbf{q})^\dagger & \mathcal{D}(\mathbf{q}) \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix}
= \varepsilon_\mu
\begin{pmatrix} X \\ Y \end{pmatrix}
}
```

where:

| Block | How to compute | Eigensystem caches needed |
|---|---|---|
| $\mathcal{A}(\mathbf{q})$ | Direct ([TDA doc §3](particle_hole.md)) | $\mathbf{k}+\mathbf{q}$ |
| $\mathcal{B}(\mathbf{q})$ | Direct (§3 of this doc) | $\mathbf{k}+\mathbf{q}$, $\mathbf{p}-\mathbf{q}$ |
| $-\mathcal{B}(\mathbf{q})^\dagger$ | From $\mathcal{B}(\mathbf{q})$ (§4.1) | no extra cache |
| $\mathcal{D}(\mathbf{q})$ | Direct (§5) or $-\mathcal{A}(-\mathbf{q})^*$ (§4.2) | $\mathbf{k}-\mathbf{q}$ |

Total eigensystem caches: **three** shifts — $+\mathbf{q}$, $-\mathbf{q}$, and $\mathbf{0}$ (the HF mesh).

> **Implementation shortcut**: When $\text{unocc}(\mathbf{k}+\mathbf{q}) = \text{unocc}(\mathbf{k}-\mathbf{q})$ for all $\mathbf{k}$ (e.g., systems with inversion symmetry, or when $\mathbf{q}$ is at a high-symmetry point), the simpler relation $\mathcal{D}(\mathbf{q}) = -\mathcal{A}(-\mathbf{q})^*$ can be used. In general, both $\mathcal{A}$ and $\mathcal{D}$ must be computed independently.
