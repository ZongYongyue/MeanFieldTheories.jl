"""
Particle-hole excitation theory (TDA / Bethe-Salpeter equation) on top of the
Hartree-Fock mean field, following the derivation in docs/src/particle_hole.md.

The central object is the effective Hamiltonian matrix (§5.7)

    H^{n0 n, n0' n'}_{k p}(q) = δ_{kp} δ_{n0 n0'} δ_{nn'} (E^n_{k+q} - E^{n0}_k)
                                + K^d_{k p}^{n0 n, n0' n'}(q)   [direct channel]
                                + K^x_{k p}^{n0 n, n0' n'}(q)   [exchange channel]

whose eigenvalues ε_μ(q) are the collective excitation energies and whose
eigenvectors ψ^{n0 n}_k(μ, q) are the particle-hole envelope functions.

Composite index: (k, j0, n)  where
  k  = excit sub-grid index (1:Nk_excit)
  j0 = index into n0_list (j0-th hole band n0 = n0_list[j0])
  n  = particle band index (from n_list, unoccupied at k+q)

Index / layout conventions (Julia is column-major):
  - eigenvectors from solve_hfk: Array{ComplexF64,3} of shape (d, d, Nk)
    evecs[:, n, ki]  ↔  U_{a,n}(k),  i.e. the n-th column at k-point ki
  - eigenvalues from solve_hfk:  Matrix{Float64} of shape (d, Nk)
    eigenvalues[n, ki]  ↔  E^n_k
"""

# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers: k-point index lookup
# ──────────────────────────────────────────────────────────────────────────────

"""
    _kindex_analytic(k, B_inv, Ns) -> Int

O(1) k-point index lookup for a **uniform rectangular** grid built by `build_kpoints`
with `box_size = (N₁, N₂, ...)` (each dimension may differ).

The grid is:  k_{m₁,m₂,...} = Σᵢ (mᵢ/Nᵢ) bᵢ,  mᵢ ∈ 0:Nᵢ-1.
`build_kpoints` iterates via `Iterators.product(0:N₁-1, 0:N₂-1, ...)`, so dim-1
varies fastest and the 1-based linear index is
  `m[1] + m[2]*N₁ + m[3]*N₁*N₂ + ... + 1`.

BZ folding is handled by `mod(round(frac[d] * Nd), Nd)`.

# Arguments
- `k`: target k-vector (may lie in any BZ image).
- `B_inv`: pinv(hcat(b₁, b₂, ...)), pre-computed from reciprocal vectors.
- `Ns::AbstractVector{Int}`: grid size along each dimension.
"""
function _kindex_analytic(k::AbstractVector, B_inv::AbstractMatrix, Ns::AbstractVector{Int})
    frac = B_inv * k                                          # fractional BZ coordinates
    m    = [mod(round(Int, frac[d] * Ns[d]), Ns[d]) for d in 1:length(Ns)]
    # linear index: dim-1 fastest (Iterators.product convention)
    idx    = m[1] + 1
    stride = Ns[1]
    for d in 2:length(Ns)
        idx    += m[d] * stride
        stride *= Ns[d]
    end
    return idx
end

# ──────────────────────────────────────────────────────────────────────────────
# Core: effective Hamiltonian matrix for a single q, all n0 in n0_list
# ──────────────────────────────────────────────────────────────────────────────

"""
    _build_heff_ph_excitation(n0_list, n_list, occ, kpoints, eigenvalues,
                              eigenvectors, kq_indices, V_k_full; excit_to_hf)
        -> (Matrix{ComplexF64}, Vector{Tuple{Int,Int,Int}})

Build the effective Hamiltonian matrix H_eff(q) for **particle-hole** excitation
momentum `q` over all hole bands `n0_list` (particle_hole.md §5.7).

The composite row/col index is `(k, j0, n)` where:
  - k  ∈ 1:Nk_excit (excit sub-grid index)
  - j0 ∈ 1:length(n0_list), with n0 = n0_list[j0]
  - n  ∈ n_list (particle band index)

A triple `(k, j0, n)` is valid when:
  1. Hole band n0 = n0_list[j0] is occupied at k:  `occ[n0, excit_to_hf[k]]`
  2. Particle band n is unoccupied at k+q:          `!occ[n, kq_indices[k]]`
  3. n ∈ n_list

# Kernel formulas (§5.6, corrected — each channel has two V terms):
  K^d = -(1/N) Σ_{abcd} [
      U*_{a,n}(k+q) U_{b,n'}(p+q) U*_{c,n0'}(p) U_{d,n0}(k)  Ṽ^{abcd}(k+q, p+q, p)
    + U*_{a,n0'}(p) U_{b,n0}(k) U*_{c,n}(k+q) U_{d,n'}(p+q)  Ṽ^{abcd}(p, k, k+q)
  ]

  K^x = +(1/N) Σ_{abcd} [
      U*_{a,n}(k+q) U_{b,n0}(k) U*_{c,n0'}(p) U_{d,n'}(p+q)  Ṽ^{abcd}(k+q, k, p)
    + U*_{a,n0'}(p) U_{b,n'}(p+q) U*_{c,n}(k+q) U_{d,n0}(k)  Ṽ^{abcd}(p, p+q, k+q)
  ]

  The two terms per channel arise from the two bilinear pairs in H_int
  (see §5.4: contractions γ-A and γ-B from Terms (I-a) and (II-a)).
  They involve V at genuinely different momentum arguments and cannot
  be combined by index relabeling when V is momentum-dependent.

# Returns
- `H::Matrix{ComplexF64}` of shape (M, M) where M = number of valid triples.
- `triples::Vector{Tuple{Int,Int,Int}}`: the `(k, j0, n)` triples.
"""
function _build_heff_ph_excitation(
    n0_list::Vector{Int},
    n_list::Vector{Int},
    occ::Matrix{Bool},
    kpoints::Vector{Vector{Float64}},
    eigenvalues::Matrix{Float64},
    eigenvectors::Array{ComplexF64, 3},
    kq_indices::Vector{Int},
    V_k_full;
    excit_to_hf::AbstractVector{Int} = 1:size(eigenvectors, 3)
)
    Nk_excit = length(excit_to_hf)
    Nn0      = length(n0_list)
    norb     = size(eigenvectors, 1)

    # ── Valid (k, j0, n) triples ───────────────────────────────────────────────
    # Note: valid n at k is determined by !occ[n, kq_indices[k]] — independent of j0.
    # Valid j0 at k requires occ[n0_list[j0], excit_to_hf[k]].
    triples = [(k, j0, n)
               for k  in 1:Nk_excit
               for j0 in 1:Nn0
               for n  in n_list
               if occ[n0_list[j0], excit_to_hf[k]] && !occ[n, kq_indices[k]]]
    M = length(triples)
    M == 0 && error("_build_heff_ph_excitation: no valid particle-hole triples for n0_list=$n0_list")

    H = zeros(ComplexF64, M, M)

    # Index lookup: for composite (k, j0), which rows of H and which n bands?
    # row_idx[(k-1)*Nn0 + j0] => Vector{Int} of row indices into H
    # n_bands[(k-1)*Nn0 + j0] => Vector{Int} of corresponding n values
    row_idx = [Int[] for _ in 1:(Nk_excit * Nn0)]
    n_bands = [Int[] for _ in 1:(Nk_excit * Nn0)]
    for (i, (k, j0, n)) in enumerate(triples)
        slot = (k - 1) * Nn0 + j0
        push!(row_idx[slot], i)
        push!(n_bands[slot], n)
    end

    # ── Diagonal: E^n_{k+q} - E^{n0}_k ───────────────────────────────────────
    for (i, (k, j0, n)) in enumerate(triples)
        n0 = n0_list[j0]
        H[i, i] = eigenvalues[n, kq_indices[k]] - eigenvalues[n0, excit_to_hf[k]]
    end

    # ── K^d + K^x: loop over (k, p) pairs, then (j0, j0') pairs ──────────────
    #
    # §5.6 kernel formulas — each channel receives TWO V contributions from
    # the two bilinear pairs of H_int (γ-A and γ-B contractions in §5.4).
    #
    # For fixed (k, p), four V tensors are computed once:
    #   Vd1 = Ṽ(k+q, p+q, p)       direct  term 1  [from (I-a)-γA]
    #   Vd2 = Ṽ(p,   k,   k+q)     direct  term 2  [from (II-a)-γB]
    #   Vx1 = Ṽ(k+q, k,   p)       exchange term 1  [from (II-a)-γA]
    #   Vx2 = Ṽ(p,   p+q, k+q)     exchange term 2  [from (I-a)-γB]
    #
    # ── Direct channel ──
    # K^d = -(1/N) [term1 + term2], where:
    #
    #   term1: Σ_{abcd} U*_an(k+q) U_bn'(p+q) U*_cn0'(p) U_dn0(k) Vd1[a,b,c,d]
    #     Free indices (a,b) → projected by U†(k+q) and U(p+q)
    #     Contracted (c,d): A_d1[a,b] = Σ_{cd} Vd1[a,b,c,d] conj(U_cn0'(p)) U_dn0(k)
    #     Vd1_r = reshape(Vd1, norb², norb²)  → rows=(a,b), cols=(c,d)
    #     kron_A = kron(U_n0_k, conj(U_n0p_p))  → index (c,d): U_n0_k[d]·conj(U_n0p_p[c])
    #
    #   term2: Σ_{abcd} U*_an0'(p) U_bn0(k) U*_cn(k+q) U_dn'(p+q) Vd2[a,b,c,d]
    #     Free indices (c,d) → projected by U†(k+q) and U(p+q)
    #     Contracted (a,b): A_d2[c,d] = Σ_{ab} Vd2[a,b,c,d] conj(U_an0'(p)) U_bn0(k)
    #     Vd2_r = reshape(permutedims(Vd2,(3,4,1,2)), norb², norb²)  → rows=(c,d), cols=(a,b)
    #     same kron_A but now contracted over (a,b): U_n0_k[b]·conj(U_n0p_p[a]) ✓
    #
    #   Both terms use the same kron vector, so: A_d = reshape((Vd1_r + Vd2_r) * kron_A, …)
    #
    # ── Exchange channel ──
    # K^x = +(1/N) [term1 + term2], where:
    #
    #   term1: Σ_{abcd} U*_an(k+q) U_bn0(k) U*_cn0'(p) U_dn'(p+q) Vx1[a,b,c,d]
    #     Free (a,d) → projected by U†(k+q) and U(p+q)
    #     Contracted (b,c): A_x1[a,d] = Σ_{bc} Vx1[a,b,c,d] U_bn0(k) conj(U_cn0'(p))
    #     Vx1_r = reshape(permutedims(Vx1,(1,4,2,3)), norb², norb²)  → rows=(a,d), cols=(b,c)
    #     kron_B = kron(conj(U_n0p_p), U_n0_k)  → index (b,c): conj(U_n0p_p[c])·U_n0_k[b]
    #
    #   term2: Σ_{abcd} U*_an0'(p) U_bn'(p+q) U*_cn(k+q) U_dn0(k) Vx2[a,b,c,d]
    #     Free (c,b) → projected by U†(k+q) [on c] and U(p+q) [on b]
    #     Contracted (a,d): A_x2[c,b] = Σ_{ad} Vx2[a,b,c,d] conj(U_an0'(p)) U_dn0(k)
    #     Vx2_r = reshape(permutedims(Vx2,(3,2,1,4)), norb², norb²)  → rows=(c,b), cols=(a,d)
    #     kron_A = kron(U_n0_k, conj(U_n0p_p))  → index (a,d): U_n0_k[d]·conj(U_n0p_p[a])
    #
    #   The two exchange terms use different kron vectors, so they are computed
    #   separately and added: A_x = A_x1 + A_x2.
    #
    # ── Final assembly ──
    # H[rows_k, rows_p] += (1/Nk) U_kq' * (A_x - A_d) * U_pq

    norb2 = norb^2

    for k in 1:Nk_excit
        kq   = kq_indices[k]
        k_hf = excit_to_hf[k]

        for p in 1:Nk_excit
            pq   = kq_indices[p]
            p_hf = excit_to_hf[p]

            # Four V tensor evaluations (one per γ contraction, §5.4)
            Vd1_raw = V_k_full(kpoints[kq], kpoints[pq],   kpoints[p_hf])  # (k+q, p+q, p)
            Vd2_raw = V_k_full(kpoints[p_hf], kpoints[k_hf], kpoints[kq]) # (p,   k,   k+q)
            Vx1_raw = V_k_full(kpoints[kq], kpoints[k_hf], kpoints[p_hf]) # (k+q, k,   p)
            Vx2_raw = V_k_full(kpoints[p_hf], kpoints[pq],  kpoints[kq])  # (p,   p+q, k+q)

            # Reshape V tensors for matrix multiplication:
            #   direct:   free indices (a,b) or (c,d), contracted = complement
            #   exchange: free indices (a,d) or (c,b), contracted = complement
            Vd1_r = reshape(Vd1_raw, norb2, norb2)                                # (a,b)×(c,d)
            Vd2_r = reshape(permutedims(Vd2_raw, (3, 4, 1, 2)), norb2, norb2)     # (c,d)×(a,b)
            Vx1_r = reshape(permutedims(Vx1_raw, (1, 4, 2, 3)), norb2, norb2)     # (a,d)×(b,c)
            Vx2_r = reshape(permutedims(Vx2_raw, (3, 2, 1, 4)), norb2, norb2)     # (c,b)×(a,d)

            # Combined direct reshaped matrix (both terms use the same kron vector)
            Vd_combined = Vd1_r + Vd2_r

            for j0 in 1:Nn0
                slot_k  = (k - 1) * Nn0 + j0
                rows_k  = row_idx[slot_k]
                isempty(rows_k) && continue

                n0      = n0_list[j0]
                U_n0_k  = eigenvectors[:, n0, k_hf]
                ns_k    = n_bands[slot_k]
                U_kq_mat = eigenvectors[:, ns_k, kq]   # norb × |ns_k|

                for j0p in 1:Nn0
                    slot_p  = (p - 1) * Nn0 + j0p
                    rows_p  = row_idx[slot_p]
                    isempty(rows_p) && continue

                    n0p      = n0_list[j0p]
                    U_n0p_p  = eigenvectors[:, n0p, p_hf]
                    ns_p     = n_bands[slot_p]
                    U_pq_mat = eigenvectors[:, ns_p, pq]  # norb × |ns_p|

                    # kron vectors for orbital-index contraction
                    # kron_A[(c/a) + (d/b - 1)*norb] = U_n0_k[d/b] · conj(U_n0p_p[c/a])
                    # kron_B[(b) + (c - 1)*norb]     = conj(U_n0p_p[c]) · U_n0_k[b]
                    kron_A = kron(U_n0_k,        conj(U_n0p_p))
                    kron_B = kron(conj(U_n0p_p), U_n0_k)

                    # Direct: A_d[free1, free2] where free = (a,b) from term1, (c,d) from term2
                    # Both project onto particle bands via U_kq' * ... * U_pq
                    A_d = reshape(Vd_combined * kron_A, norb, norb)

                    # Exchange: two terms with different kron vectors
                    A_x1 = reshape(Vx1_r * kron_B, norb, norb)   # free=(a,d)
                    A_x2 = reshape(Vx2_r * kron_A, norb, norb)   # free=(c,b)
                    A_x  = A_x1 + A_x2

                    H[rows_k, rows_p] .+= (U_kq_mat' * (A_x - A_d) * U_pq_mat) ./ Nk_excit
                end
            end
        end
    end

    return H, triples
end

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

"""
    solve_ph_excitations(dofs, twobody, hf_result, qpoints, reciprocal_vecs;
                         excit_kpoints=nothing, n0_list=nothing, n_list=nothing,
                         verbose=true) -> NamedTuple

Compute the **particle-hole** excitation spectrum (TDA / Bethe-Salpeter) on top
of a converged Hartree-Fock ground state.

For each excitation momentum q ∈ `qpoints`, builds the unified effective
Hamiltonian H_eff(q) with composite index `(k, j0, n)` over all hole bands
`n0_list` simultaneously, then diagonalizes it.

# Arguments
- `dofs::SystemDofs`: internal DOFs of one magnetic unit cell.
- `twobody`: two-body interaction terms; needs `.ops` and `.irvec`.
- `hf_result`: NamedTuple returned by `solve_hfk`.  Fields used:
  `eigenvalues` (d×Nk_hf), `eigenvectors` (d×d×Nk_hf), `kpoints`, `mu`.
- `qpoints::Vector{Vector{Float64}}`: excitation momenta to compute.
- `reciprocal_vecs::Vector{Vector{Float64}}`: reciprocal lattice vectors b_i.

# Keyword Arguments
- `excit_kpoints`: coarser k-grid for the particle-hole sum.  `nothing` uses
  the full HF k-grid.
- `n0_list`: hole band indices.  `nothing` = all bands occupied at any excit k.
- `n_list`: particle band indices.  `nothing` = all bands (occupancy filter applied).
- `verbose::Bool = true`: print progress information.

# Returns
NamedTuple with fields:
- `qpoints`: the input q-points (echoed).
- `n0_list::Vector{Int}`: the hole bands actually used.
- `n_list::Vector{Int}`: the particle bands actually used.
- `energies::Vector{Vector{Float64}}`: `energies[qi]` is the sorted vector of
  eigenvalues ε_μ(q) for q-point index qi (all n0 channels unified).
- `wavefunctions::Vector{Matrix{ComplexF64}}`: `wavefunctions[qi]` is the
  eigenvector matrix; column μ is the envelope function ψ(μ) with composite
  index `(k, j0, n)`.
- `triples::Vector{Vector{Tuple{Int,Int,Int}}}`: `triples[qi]` lists the
  `(k, j0, n)` triples for each row/column of `wavefunctions[qi]`.
"""
function solve_ph_excitations(
    dofs::SystemDofs,
    twobody,
    hf_result::NamedTuple,
    qpoints::Vector{Vector{Float64}},
    reciprocal_vecs::Vector{Vector{Float64}};
    excit_kpoints::Union{Nothing, Vector{Vector{Float64}}} = nothing,
    n0_list::Union{Nothing, Vector{Int}} = nothing,
    n_list::Union{Nothing, Vector{Int}}  = nothing,
    verbose::Bool = true
)
    eigenvalues  = hf_result.eigenvalues   # (norb, Nk_hf)
    eigenvectors = hf_result.eigenvectors  # (norb, norb, Nk_hf)
    kpoints_hf   = hf_result.kpoints       # Vector{Vector{Float64}}, length Nk_hf
    mu           = hf_result.mu

    norb    = size(eigenvalues, 1)
    tol_occ = 1e-8

    # ── Detect active dimensions in HF k-point grid ───────────────────────────
    D_hf  = length(kpoints_hf[1])
    Nk_hf = length(kpoints_hf)
    active_dims = [d for d in 1:D_hf if any(abs(k[d]) > 1e-12 for k in kpoints_hf)]
    n_active    = length(active_dims)

    length(reciprocal_vecs) == n_active ||
        error("reciprocal_vecs must have $n_active vector(s) for this k-point grid " *
              "(detected $n_active active dimension(s) in hf.kpoints), " *
              "got $(length(reciprocal_vecs))")

    # active_k: project any k-vector (D_hf-dim or n_active-dim) onto active subspace
    function active_k(k::AbstractVector)
        length(k) == n_active && return collect(Float64, k)
        length(k) == D_hf     && return Float64[k[d] for d in active_dims]
        error("k length $(length(k)) incompatible: expected $n_active or $D_hf")
    end

    # B_inv in active subspace (n_active × n_active)
    rv_active = [active_k(b) for b in reciprocal_vecs]
    B_inv     = inv(hcat(rv_active...))

    # Infer per-dimension grid sizes from fractional coordinates in active subspace
    Ns_hf = let frac_sets = [Set{Int}() for _ in 1:n_active]
        for k in kpoints_hf
            f = B_inv * active_k(k)
            for d in 1:n_active
                push!(frac_sets[d], round(Int, f[d] * Nk_hf))
            end
        end
        [length(s) for s in frac_sets]
    end
    @assert prod(Ns_hf) == Nk_hf "HF k-grid size mismatch: inferred $(Ns_hf) but got $Nk_hf points"

    # ── Excit sub-grid (defaults to full HF grid) ─────────────────────────────
    if excit_kpoints === nothing
        kpoints_excit = kpoints_hf
        excit_to_hf   = collect(1:Nk_hf)
    else
        kpoints_excit = excit_kpoints
        excit_to_hf   = [_kindex_analytic(active_k(k), B_inv, Ns_hf) for k in kpoints_excit]
    end
    Nk_excit = length(kpoints_excit)

    # ── Occupied / unoccupied band masks (on full HF grid) ───────────────────
    occ = Matrix{Bool}(eigenvalues .<= mu + tol_occ)   # (norb, Nk_hf)

    # ── Resolve hole bands: occupied at any excit k-point ────────────────────
    if n0_list === nothing
        n0_list = [n for n in 1:norb if any(occ[n, excit_to_hf])]
    end
    isempty(n0_list) && error("solve_ph_excitations: no occupied bands found")

    # ── Resolve particle bands ────────────────────────────────────────────────
    if n_list === nothing
        n_list = collect(1:norb)
    end

    # ── Build full interaction kernel ─────────────────────────────────────────
    V_k_full = build_Vk(build_Vr(dofs, twobody.ops, twobody.irvec))

    # ── Output arrays ─────────────────────────────────────────────────────────
    Nq            = length(qpoints)
    energies      = Vector{Vector{Float64}}(undef, Nq)
    wavefunctions = Vector{Matrix{ComplexF64}}(undef, Nq)
    triples_out   = Vector{Vector{Tuple{Int,Int,Int}}}(undef, Nq)

    # ── Main loop: parallel over q-points ────────────────────────────────────
    print_lock = verbose ? ReentrantLock() : nothing
    Threads.@threads for qi in 1:Nq
        q = qpoints[qi]
        if verbose
            lock(print_lock) do
                println("q-point $qi / $Nq : q = $q")
            end
        end

        # kq_indices[k]: HF index of (excit_k + q), projected onto active subspace
        kq_indices = [_kindex_analytic(active_k(kpoints_excit[k]) .+ active_k(q), B_inv, Ns_hf)
                      for k in 1:Nk_excit]

        H, triples = _build_heff_ph_excitation(
            n0_list, n_list, occ, kpoints_hf, eigenvalues, eigenvectors,
            kq_indices, V_k_full; excit_to_hf = excit_to_hf
        )
        F = eigen(Hermitian(H))
        energies[qi]      = F.values
        wavefunctions[qi] = F.vectors
        triples_out[qi]   = triples
    end

    return (
        qpoints       = qpoints,
        n0_list       = n0_list,
        n_list        = n_list,
        energies      = energies,
        wavefunctions = wavefunctions,
        triples       = triples_out,
    )
end