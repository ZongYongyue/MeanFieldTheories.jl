"""
Particle-hole excitation theory (TDA / RPA) on top of the Hartree-Fock mean
field, following the derivation in docs/src/particle_hole.md and
docs/src/particle_hole2.md.

## TDA (solver=:TDA, default)
Diagonalizes the M×M Hermitian effective Hamiltonian A(q).

## RPA (solver=:RPA)
Solves the full 2M×2M non-Hermitian eigenvalue problem:
    (A  B; -B†  D)(X; Y) = ε(X; Y)
where C = -B† (§4.1 of particle_hole2.md) and D is the backward-backward
block.  When ±q share the same occupation spaces, D(q) = -A(-q)* is used
(§2.6); otherwise D is computed directly (§5).

## Momentum handling
Shifted momenta (k+q, p-q, …) are first folded back into the first BZ.
If the folded point coincides with a k-mesh point, its HF eigensystem is reused
directly.  Otherwise, the HF Hamiltonian is rebuilt and diagonalized on the fly
via `energy_bands`.

## Index conventions
  - eigenvectors: (d, d, Nk);  evecs[:, n, ki] ↔ U_{a,n}(k)
  - eigenvalues:  (d, Nk);     eigenvalues[n, ki] ↔ E^n_k
"""


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers — BZ folding & eigensystem cache
# ──────────────────────────────────────────────────────────────────────────────

"""
    _fold_to_BZ(k, B_inv, B_mat) -> Vector{Float64}

Fold an arbitrary k-vector back into the first Brillouin zone by mapping its
fractional coordinates to [0, 1).
"""
function _fold_to_BZ(k::AbstractVector, B_inv::AbstractMatrix, B_mat::AbstractMatrix)
    frac = B_inv * k
    frac_folded = [mod(f, 1.0) for f in frac]
    # Snap values very close to 1.0 back to 0.0
    for i in eachindex(frac_folded)
        if frac_folded[i] > 1.0 - 1e-10
            frac_folded[i] = 0.0
        end
    end
    return B_mat * frac_folded
end

"""
    _find_on_mesh(k_folded, kpoints_folded; tol) -> Union{Int, Nothing}

Check if `k_folded` matches any point in `kpoints_folded` within tolerance.
Returns the index if found, `nothing` otherwise.
"""
function _find_on_mesh(k_folded::AbstractVector, kpoints_folded::Vector{Vector{Float64}};
                       tol::Real = 1e-8)
    for (i, kp) in enumerate(kpoints_folded)
        if all(abs.(k_folded .- kp) .< tol)
            return i
        end
    end
    return nothing
end

"""
    EigenCache

Caches eigenvalues and eigenvectors at arbitrary k-points.
On-mesh points reuse the HF result; off-mesh points are computed via `energy_bands`.
"""
struct EigenCache
    evals::Matrix{Float64}          # (norb, N_total) — all cached eigenvalues
    evecs::Array{ComplexF64, 3}     # (norb, norb, N_total) — all cached eigenvectors
    kvecs::Vector{Vector{Float64}}  # (N_total,) — the actual (folded) k-vectors
    index_map::Dict{Int, Int}       # mesh ki -> cache index (for on-mesh lookups)
end

"""
    _build_eigen_cache(shift, kpoints_hf, kpoints_folded, B_inv, B_mat,
                       evals_hf, evecs_hf, dofs, onebody, twobody, G_k)

For each k in `kpoints_hf`, compute eigensystem at `k + shift`:
  1. Fold `k + shift` into the first BZ.
  2. If the folded point is on the HF mesh, reuse its eigensystem.
  3. Otherwise, batch-compute via `energy_bands`.

Returns an EigenCache with Nk entries (one per HF k-point).
"""
function _build_eigen_cache(
    shift::AbstractVector,
    kpoints_hf::Vector{Vector{Float64}},
    kpoints_folded::Vector{Vector{Float64}},
    B_inv::AbstractMatrix,
    B_mat::AbstractMatrix,
    evals_hf::Matrix{Float64},
    evecs_hf::Array{ComplexF64, 3},
    dofs::SystemDofs,
    onebody,
    twobody,
    G_k::Array{ComplexF64, 3}
)
    Nk   = length(kpoints_hf)
    norb = size(evals_hf, 1)

    # For each ki, fold k_i + shift and check if it's on mesh
    folded_targets = Vector{Vector{Float64}}(undef, Nk)
    mesh_hits      = Vector{Union{Int, Nothing}}(undef, Nk)  # index into HF mesh, or nothing
    off_mesh_ki    = Int[]                                     # ki indices needing computation

    for ki in 1:Nk
        k_shifted = kpoints_hf[ki] .+ shift
        k_fold    = _fold_to_BZ(k_shifted, B_inv, B_mat)
        folded_targets[ki] = k_fold
        hit = _find_on_mesh(k_fold, kpoints_folded)
        mesh_hits[ki] = hit
        if hit === nothing
            push!(off_mesh_ki, ki)
        end
    end

    # Initialize cache arrays
    cache_evals = Matrix{Float64}(undef, norb, Nk)
    cache_evecs = Array{ComplexF64, 3}(undef, norb, norb, Nk)

    # Fill on-mesh entries from HF data
    index_map = Dict{Int, Int}()
    for ki in 1:Nk
        hit = mesh_hits[ki]
        if hit !== nothing
            cache_evals[:, ki]    = evals_hf[:, hit]
            cache_evecs[:, :, ki] = evecs_hf[:, :, hit]
        end
        index_map[ki] = ki
    end

    # Batch-compute off-mesh entries
    if !isempty(off_mesh_ki)
        off_mesh_kpts = [folded_targets[ki] for ki in off_mesh_ki]
        bands_off, evecs_off = energy_bands(dofs, onebody, twobody,
                                            kpoints_hf, G_k, off_mesh_kpts)
        for (j, ki) in enumerate(off_mesh_ki)
            cache_evals[:, ki]    = bands_off[:, j]
            cache_evecs[:, :, ki] = evecs_off[:, :, j]
        end
    end

    return EigenCache(cache_evals, cache_evecs,
                      folded_targets, index_map)
end


# ──────────────────────────────────────────────────────────────────────────────
# Particle-hole pair construction
# ──────────────────────────────────────────────────────────────────────────────

"""
    _build_ph_pairs(evals_k, evals_kq, mu; tol_occ, n_list) -> Vector{Tuple{Int,Int,Int}}

Build the list of valid particle-hole triples (ki, n0, n) for a given q.
Occupation is determined automatically per k-point:
  - n0 ∈ occ(k):     E^{n0}_k < mu
  - n  ∈ unocc(k+q): E^n_{k+q} > mu  AND  n ∈ n_list (if provided)
"""
function _build_ph_pairs(
    evals_k::AbstractMatrix,
    evals_kq::AbstractMatrix,
    mu::Real;
    tol_occ::Real = 1e-8,
    n_list::Union{Nothing, Vector{Int}} = nothing
)
    norb = size(evals_k, 1)
    Nk   = size(evals_k, 2)
    candidates = n_list === nothing ? (1:norb) : n_list

    triples = Tuple{Int,Int,Int}[]
    for ki in 1:Nk
        occ   = [n for n in 1:norb    if evals_k[n, ki]  < mu - tol_occ]
        unocc = [n for n in candidates if evals_kq[n, ki] > mu + tol_occ]
        for n0 in occ, n in unocc
            push!(triples, (ki, n0, n))
        end
    end
    return triples
end


# ──────────────────────────────────────────────────────────────────────────────
# A matrix (TDA effective Hamiltonian)
# ──────────────────────────────────────────────────────────────────────────────

"""
    _build_A_matrix(V_k, evecs_k, evals_k, cache_kq,
                    kpoints, q, triples, Nk) -> Matrix{ComplexF64}

Build the TDA effective Hamiltonian A(q), following §3 of particle_hole.md:

    A^{n₀n, n₀'n'}_{kp}(q) = δ_{kp} δ_{n₀n₀'} δ_{nn'} (E^n_{k+q} - E^{n₀}_k)
                              + Kd + Kx

Uses kron+reshape+BLAS for fast tensor contractions.
"""
function _build_A_matrix(
    V_k,
    evecs_k::AbstractArray{<:Number,3},
    evals_k::AbstractMatrix,
    cache_kq::EigenCache,
    kpoints::Vector{Vector{Float64}},
    q::Vector{Float64},
    triples::Vector{Tuple{Int,Int,Int}},
    Nk::Int
)
    M    = length(triples)
    norb = size(evals_k, 1)
    norb2 = norb * norb
    A    = zeros(ComplexF64, M, M)

    # Diagonal: free particle-hole pair energy
    for (I, (ki, n0, n)) in enumerate(triples)
        A[I, I] = cache_kq.evals[n, ki] - evals_k[n0, ki]
    end

    # ── Build index lookup: group triples by (ki, n0) ──
    # slot_key[ki][n0] -> (row_indices, n_bands)
    slot_key = Dict{Tuple{Int,Int}, Tuple{Vector{Int}, Vector{Int}}}()
    for (I, (ki, n0, n)) in enumerate(triples)
        key = (ki, n0)
        if !haskey(slot_key, key)
            slot_key[key] = (Int[], Int[])
        end
        push!(slot_key[key][1], I)
        push!(slot_key[key][2], n)
    end

    # Collect unique ki and n0 values per ki
    ki_set = sort(unique(t[1] for t in triples))
    n0_per_ki = Dict{Int, Vector{Int}}()
    for ki in ki_set
        n0_per_ki[ki] = sort(unique(t[2] for t in triples if t[1] == ki))
    end

    # ── Interaction kernel ──
    invN = 1.0 / Nk
    for ki in ki_set
        kk  = kpoints[ki]
        kkq = kk .+ q

        Uk  = @view evecs_k[:, :, ki]
        Ukq = @view cache_kq.evecs[:, :, ki]

        for pi in ki_set
            kp  = kpoints[pi]
            kpq = kp .+ q

            Up  = @view evecs_k[:, :, pi]
            Upq = @view cache_kq.evecs[:, :, pi]

            # ── Four V calls ──
            Vd1_raw = V_k(kkq, kpq, kp)     # (k+q, p+q, p)
            Vd2_raw = V_k(kp, kk, kkq)      # (p,   k,   k+q)
            Vx1_raw = V_k(kkq, kk, kp)      # (k+q, k,   p)
            Vx2_raw = V_k(kp, kpq, kkq)     # (p,   p+q, k+q)

            # ── Reshape V for A kernel ──
            Vd_A  = reshape(Vd1_raw, norb2, norb2) .+
                     reshape(permutedims(Vd2_raw, (3,4,1,2)), norb2, norb2)
            Vx1_A = reshape(permutedims(Vx1_raw, (1,4,2,3)), norb2, norb2)
            Vx2_A = reshape(permutedims(Vx2_raw, (3,2,1,4)), norb2, norb2)

            # ── Inner loops over hole bands ──
            for n0 in get(n0_per_ki, ki, Int[])
                key_k = (ki, n0)
                haskey(slot_key, key_k) || continue
                rows_k, ns_k = slot_key[key_k]

                U_n0_k   = @view Uk[:, n0]
                U_kq_mat = @view Ukq[:, ns_k]

                for n0p in get(n0_per_ki, pi, Int[])
                    key_p = (pi, n0p)
                    haskey(slot_key, key_p) || continue
                    rows_p, ns_p = slot_key[key_p]

                    U_n0p_p  = @view Up[:, n0p]
                    U_pq_mat = @view Upq[:, ns_p]

                    conj_U_n0p_p = conj(U_n0p_p)
                    kron_Ac = kron(U_n0_k,       conj_U_n0p_p)
                    kron_Ax = kron(conj_U_n0p_p, U_n0_k)

                    Ad = reshape(Vd_A  * kron_Ac, norb, norb)
                    Ax = reshape(Vx1_A * kron_Ax, norb, norb) +
                         reshape(Vx2_A * kron_Ac, norb, norb)

                    A[rows_k, rows_p] .+= (U_kq_mat' * (Ax - Ad) * U_pq_mat) .* invN
                end
            end
        end
    end

    return A
end


# ──────────────────────────────────────────────────────────────────────────────
# B matrix (RPA coupling)
# ──────────────────────────────────────────────────────────────────────────────

"""
    _build_B_matrix(V_k, evecs_k, cache_kq, cache_kmq,
                    kpoints, q, triples, Nk) -> Matrix{ComplexF64}

Build the RPA coupling matrix B(q), following §3.5 of particle_hole2.md:

    B^{n₀n, n₀'n'}_{kp}(q) = Bd + Bx

Uses eigensystems at k+q (from `cache_kq`) and p-q (from `cache_kmq`).
Uses kron+reshape+BLAS for fast tensor contractions.
"""
function _build_B_matrix(
    V_k,
    evecs_k::AbstractArray{<:Number,3},
    cache_kq::EigenCache,
    cache_kmq::EigenCache,
    kpoints::Vector{Vector{Float64}},
    q::Vector{Float64},
    triples::Vector{Tuple{Int,Int,Int}},
    Nk::Int
)
    M    = length(triples)
    norb = size(evecs_k, 1)
    norb2 = norb * norb
    B    = zeros(ComplexF64, M, M)

    # ── Build index lookup: group triples by (ki, n0) ──
    slot_key = Dict{Tuple{Int,Int}, Tuple{Vector{Int}, Vector{Int}}}()
    for (I, (ki, n0, n)) in enumerate(triples)
        key = (ki, n0)
        if !haskey(slot_key, key)
            slot_key[key] = (Int[], Int[])
        end
        push!(slot_key[key][1], I)
        push!(slot_key[key][2], n)
    end

    ki_set = sort(unique(t[1] for t in triples))
    n0_per_ki = Dict{Int, Vector{Int}}()
    for ki in ki_set
        n0_per_ki[ki] = sort(unique(t[2] for t in triples if t[1] == ki))
    end

    # ── Interaction kernel ──
    invN = 1.0 / Nk
    for ki in ki_set
        kk  = kpoints[ki]
        kkq = kk .+ q

        Uk   = @view evecs_k[:, :, ki]
        Ukq  = @view cache_kq.evecs[:, :, ki]

        for pi in ki_set
            kp   = kpoints[pi]
            kpmq = kp .- q

            Up   = @view evecs_k[:, :, pi]
            Upmq = @view cache_kmq.evecs[:, :, pi]

            # ── Four V calls for B ──
            Vd1_raw = V_k(kpmq, kk, kkq)     # (p-q, k,   k+q)
            Vd2_raw = V_k(kkq, kp, kpmq)     # (k+q, p,   p-q)
            Vx1_raw = V_k(kpmq, kp, kkq)     # (p-q, p,   k+q)
            Vx2_raw = V_k(kkq, kk, kpmq)     # (k+q, k,   p-q)

            # ── Reshape V for B kernel ──
            # Bd1: V[a,b,c,d] with n↔c, n'↔a, n0↔b, n0'↔d → perm (3,1,2,4)
            # Bd2: V[a,b,c,d] with n↔a, n'↔c, n0'↔b, n0↔d → perm (1,3,2,4)
            # Bx1: V[a,b,c,d] with n↔c, n'↔a, n0'↔b, n0↔d → perm (3,1,2,4)
            # Bx2: V[a,b,c,d] with n↔a, n'↔c, n0↔b, n0'↔d → perm (1,3,2,4)
            Vd1_B = reshape(permutedims(Vd1_raw, (3,1,2,4)), norb2, norb2)
            Vd2_B = reshape(permutedims(Vd2_raw, (1,3,2,4)), norb2, norb2)
            Vx1_B = reshape(permutedims(Vx1_raw, (3,1,2,4)), norb2, norb2)
            Vx2_B = reshape(permutedims(Vx2_raw, (1,3,2,4)), norb2, norb2)

            # ── Inner loops over hole bands ──
            for n0 in get(n0_per_ki, ki, Int[])
                key_k = (ki, n0)
                haskey(slot_key, key_k) || continue
                rows_k, ns_k = slot_key[key_k]

                U_n0_k   = @view Uk[:, n0]
                U_kq_mat = @view Ukq[:, ns_k]

                for n0p in get(n0_per_ki, pi, Int[])
                    key_p = (pi, n0p)
                    haskey(slot_key, key_p) || continue
                    rows_p, ns_p = slot_key[key_p]

                    U_n0p_p   = @view Up[:, n0p]
                    U_pmq_mat = @view Upmq[:, ns_p]

                    # kron vectors for hole contractions
                    # Bd1, Bx2: n0↔b, n0'↔d → kron(U_n0_k, U_n0p_p)
                    # Bd2, Bx1: n0'↔b, n0↔d → kron(U_n0p_p, U_n0_k)
                    kron_B1 = kron(U_n0_k,  U_n0p_p)
                    kron_B2 = kron(U_n0p_p, U_n0_k)

                    Bd = reshape(Vd1_B * kron_B1, norb, norb) +
                         reshape(Vd2_B * kron_B2, norb, norb)
                    Bx = reshape(Vx1_B * kron_B2, norb, norb) +
                         reshape(Vx2_B * kron_B1, norb, norb)

                    # Project: left = U*(k+q), right = U*(p-q)
                    B[rows_k, rows_p] .+= (U_kq_mat' * (Bx - Bd) * conj(U_pmq_mat)) .* invN
                end
            end
        end
    end

    return B
end


# ──────────────────────────────────────────────────────────────────────────────
# D matrix (RPA backward-backward block)
# ──────────────────────────────────────────────────────────────────────────────

"""
    _build_D_matrix(V_k, evecs_k, evals_k, cache_kmq,
                    kpoints, q, triples, mu, Nk) -> Matrix{ComplexF64}

Build the RPA backward-backward matrix D(q), following §5.5–5.6 of
particle_hole2.md:

    D^{n₀n, n₀'n'}_{kp}(q) = δ_{kp} δ_{n₀n₀'} δ_{nn'} (E^{n₀}_k - E^n_{k-q})
                              + θ_n(k-q) θ_{n'}(p-q) × (Dd + Dx)

where θ_n(k-q) = 1 if n ∈ unocc(k-q), 0 otherwise.

Uses eigensystems at k-q (from `cache_kmq`).  Note the opposite sign structure
compared to A: direct +1/N, exchange -1/N.
Uses kron+reshape+BLAS for fast tensor contractions.
"""
function _build_D_matrix(
    V_k,
    evecs_k::AbstractArray{<:Number,3},
    evals_k::AbstractMatrix,
    cache_kmq::EigenCache,
    kpoints::Vector{Vector{Float64}},
    q::Vector{Float64},
    triples::Vector{Tuple{Int,Int,Int}},
    mu::Real,
    Nk::Int;
    tol_occ::Real = 1e-8
)
    M    = length(triples)
    norb = size(evals_k, 1)
    norb2 = norb * norb
    D    = zeros(ComplexF64, M, M)

    # Diagonal: E^{n₀}_k - E^n_{k-q}
    for (I, (ki, n0, n)) in enumerate(triples)
        D[I, I] = evals_k[n0, ki] - cache_kmq.evals[n, ki]
    end

    # ── Build index lookup: group triples by (ki, n0) ──
    slot_key = Dict{Tuple{Int,Int}, Tuple{Vector{Int}, Vector{Int}}}()
    for (I, (ki, n0, n)) in enumerate(triples)
        key = (ki, n0)
        if !haskey(slot_key, key)
            slot_key[key] = (Int[], Int[])
        end
        push!(slot_key[key][1], I)
        push!(slot_key[key][2], n)
    end

    ki_set = sort(unique(t[1] for t in triples))
    n0_per_ki = Dict{Int, Vector{Int}}()
    for ki in ki_set
        n0_per_ki[ki] = sort(unique(t[2] for t in triples if t[1] == ki))
    end

    # ── Interaction kernel ──
    invN = 1.0 / Nk
    for ki in ki_set
        kk   = kpoints[ki]
        kkmq = kk .- q

        Uk   = @view evecs_k[:, :, ki]
        Ukmq = @view cache_kmq.evecs[:, :, ki]

        for pi in ki_set
            kp   = kpoints[pi]
            kpmq = kp .- q

            Up   = @view evecs_k[:, :, pi]
            Upmq = @view cache_kmq.evecs[:, :, pi]

            # ── Four V calls for D ──
            Vd1_raw = V_k(kpmq, kkmq, kk)    # (p-q, k-q, k)
            Vd2_raw = V_k(kk, kp, kpmq)      # (k,   p,   p-q)
            Vx1_raw = V_k(kpmq, kp, kk)      # (p-q, p,   k)
            Vx2_raw = V_k(kk, kkmq, kpmq)    # (k,   k-q, p-q)

            # ── Reshape V for D kernel ──
            # Dd1: rows=(b,a) particle, cols=(c,d) hole → perm (2,1,3,4)
            # Dd2: rows=(d,c) particle, cols=(a,b) hole → perm (4,3,1,2)
            # Dx1: rows=(d,a) particle, cols=(b,c) hole → perm (4,1,2,3)
            # Dx2: rows=(b,c) particle, cols=(a,d) hole → perm (2,3,1,4)
            Vd1_D = reshape(permutedims(Vd1_raw, (2,1,3,4)), norb2, norb2)
            Vd2_D = reshape(permutedims(Vd2_raw, (4,3,1,2)), norb2, norb2)
            Vx1_D = reshape(permutedims(Vx1_raw, (4,1,2,3)), norb2, norb2)
            Vx2_D = reshape(permutedims(Vx2_raw, (2,3,1,4)), norb2, norb2)

            # ── Inner loops over hole bands ──
            for n0 in get(n0_per_ki, ki, Int[])
                key_k = (ki, n0)
                haskey(slot_key, key_k) || continue
                rows_k, ns_k = slot_key[key_k]

                # Occupation factor θ_n(k-q) for each n in ns_k
                theta_k = [cache_kmq.evals[n, ki] > mu + tol_occ ? 1.0 : 0.0 for n in ns_k]
                all(theta_k .== 0.0) && continue

                U_n0_k    = @view Uk[:, n0]
                U_kmq_mat = @view Ukmq[:, ns_k]

                for n0p in get(n0_per_ki, pi, Int[])
                    key_p = (pi, n0p)
                    haskey(slot_key, key_p) || continue
                    rows_p, ns_p = slot_key[key_p]

                    # Occupation factor θ_{n'}(p-q) for each n' in ns_p
                    theta_p = [cache_kmq.evals[n, pi] > mu + tol_occ ? 1.0 : 0.0 for n in ns_p]
                    all(theta_p .== 0.0) && continue

                    U_n0p_p   = @view Up[:, n0p]
                    U_pmq_mat = @view Upmq[:, ns_p]

                    # kron vectors for hole contractions
                    # Dd1, Dd2, Dx2: kron(conj(U_n0_k), U_n0p_p)
                    # Dx1:           kron(U_n0p_p, conj(U_n0_k))
                    conj_U_n0_k = conj(U_n0_k)
                    kron_Dc1 = kron(conj_U_n0_k, U_n0p_p)
                    kron_Dc2 = kron(U_n0p_p, conj_U_n0_k)

                    Dd = reshape(Vd1_D * kron_Dc1, norb, norb) +
                         reshape(Vd2_D * kron_Dc1, norb, norb)
                    Dx = reshape(Vx1_D * kron_Dc2, norb, norb) +
                         reshape(Vx2_D * kron_Dc1, norb, norb)

                    # Project: left = U*(k-q), right = U*(p-q)
                    # D has opposite sign: +(direct) -(exchange)
                    block = (U_kmq_mat' * (Dd - Dx) * conj(U_pmq_mat)) .* invN

                    # Apply occupation factors θ_n(k-q) and θ_{n'}(p-q)
                    block .*= theta_k       # multiply rows
                    block .*= theta_p'      # multiply columns

                    D[rows_k, rows_p] .+= block
                end
            end
        end
    end

    return D
end


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

"""
    solve_ph_excitations(dofs, onebody, twobody, hf_result, qpoints,
                         reciprocal_vecs; solver=:TDA, n_list=nothing,
                         eta=1e-10, verbose=true) -> NamedTuple

Compute the particle-hole excitation spectrum along a q-path.

Shifted momenta (k+q, p-q, …) are folded into the first BZ.  If the folded
point falls on the HF k-mesh, its eigensystem is reused; otherwise it is
recomputed on the fly via `energy_bands`.

Particle-hole pairs are auto-detected per k-point from eigenvalues vs mu.

# Keyword `n_list`
Optional vector of band indices allowed as particle (unoccupied) states.
Default `nothing` means all bands are candidates.  Use e.g. `n_list=[1,2]`
to restrict the particle-hole space to low-energy bands only.

# Keyword `eta`
Regularization added to the diagonal of the Bosonic BdG matrix before
Cholesky decomposition (RPA only).  Default `1e-10`.  Increase (e.g. `1e-6`)
for systems where the matrix is nearly singular near Goldstone modes.

# Keyword `solver`
- `:TDA` (default): Tamm-Dancoff (Hermitian, M×M). Builds only A.
- `:RPA`: Full RPA. Builds A, B, and D, then solves the 2M×2M eigenvalue
  problem:  (A  B; -B†  D)(X; Y) = ε(X; Y).
  The D block uses D(q) = -A(-q)* when the ±q occupation spaces match,
  otherwise it is computed directly with the residual occupation factor.
"""
function solve_ph_excitations(
    dofs::SystemDofs,
    onebody,
    twobody,
    hf_result::NamedTuple,
    qpoints::Vector{Vector{Float64}},
    reciprocal_vecs::Vector{Vector{Float64}};
    solver::Symbol = :TDA,
    n_list::Union{Nothing, Vector{Int}} = nothing,
    eta::Float64   = 1e-10,
    verbose::Bool  = true
)
    solver in (:TDA, :RPA) || error("solver must be :TDA or :RPA, got :$solver")

    evals_hf     = hf_result.eigenvalues
    evecs_hf     = hf_result.eigenvectors
    kpoints_hf   = hf_result.kpoints
    G_k          = hf_result.G_k
    mu           = hf_result.mu
    Nk           = length(kpoints_hf)
    tol_occ      = 1e-8

    # Reciprocal-space geometry (pinv handles D_lat < D_embed, e.g. 1D chain in 2D space)
    B_mat = hcat(reciprocal_vecs...)
    B_inv = pinv(B_mat)

    # Precompute folded HF k-points for mesh-hit detection
    kpoints_folded = [_fold_to_BZ(k, B_inv, B_mat) for k in kpoints_hf]

    # Build the interaction kernel in momentum space
    V_r = build_Vr(dofs, twobody.ops, twobody.irvec)
    V_k = build_Vk(V_r)

    Nq = length(qpoints)
    energies      = Vector{Vector{Float64}}(undef, Nq)
    wavefunctions = Vector{Matrix{ComplexF64}}(undef, Nq)
    triples_out   = Vector{Vector{Tuple{Int,Int,Int}}}(undef, Nq)

    verbose && println("  Solving $Nq q-points ($solver) with $(Threads.nthreads()) thread(s)...")
    log_lock = ReentrantLock()

    Threads.@threads for iq in 1:Nq
        q = qpoints[iq]

        # Build eigensystem cache for k+q shift
        cache_kq = _build_eigen_cache(q, kpoints_hf, kpoints_folded,
                                      B_inv, B_mat, evals_hf, evecs_hf,
                                      dofs, onebody, twobody, G_k)

        # Build particle-hole pairs (occupation auto-detected per k-point)
        triples = _build_ph_pairs(evals_hf, cache_kq.evals, mu;
                                  tol_occ=tol_occ, n_list=n_list)
        M = length(triples)
        triples_out[iq] = triples

        if M == 0
            energies[iq]      = Float64[]
            wavefunctions[iq] = zeros(ComplexF64, 0, 0)
            verbose && lock(log_lock) do
                @printf("  q[%d/%d] = [%s] ... M=0, skip\n",
                    iq, Nq, join([@sprintf("%.4f", qi) for qi in q], ", "))
            end
            continue
        end

        # Build the A matrix (TDA effective Hamiltonian)
        A = _build_A_matrix(V_k, evecs_hf, evals_hf, cache_kq,
                            kpoints_hf, q, triples, Nk)

        if solver == :TDA
            F = eigen(Hermitian(A))
            energies[iq]      = F.values
            wavefunctions[iq] = F.vectors
            verbose && lock(log_lock) do
                @printf("  q[%d/%d] = [%s] ... M=%d, done (TDA)\n",
                    iq, Nq, join([@sprintf("%.4f", qi) for qi in q], ", "), M)
            end

        elseif solver == :RPA
            # Build eigensystem cache for k-q shift (needed by B and D)
            mq = -q
            cache_kmq = _build_eigen_cache(mq, kpoints_hf, kpoints_folded,
                                           B_inv, B_mat, evals_hf, evecs_hf,
                                           dofs, onebody, twobody, G_k)

            # Build B(q): needs cache_kq (k+q) and cache_kmq (p-q)
            B_q = _build_B_matrix(V_k, evecs_hf, cache_kq, cache_kmq,
                                  kpoints_hf, q, triples, Nk)

            # Check if +q and -q triples match (for D shortcut)
            triples_mq = _build_ph_pairs(evals_hf, cache_kmq.evals, mu;
                                         tol_occ=tol_occ, n_list=n_list)
            same_pairs = (triples == triples_mq)
            d_method = ""

            if same_pairs
                # D(q) = -A(-q)* — build A at -q and negate-conjugate
                A_mq = _build_A_matrix(V_k, evecs_hf, evals_hf, cache_kmq,
                                       kpoints_hf, mq, triples, Nk)
                D_q = -conj.(A_mq)
                d_method = "D=-A(-q)*"
            else
                # Compute D(q) directly with occupation factor
                D_q = _build_D_matrix(V_k, evecs_hf, evals_hf, cache_kmq,
                                      kpoints_hf, q, triples, mu, Nk;
                                      tol_occ=tol_occ)
                d_method = "D=direct"
            end

            # ── Bosonic BdG Cholesky method (Shindou et al. PRB 87, 174427, 2013) ──
            cholesky_ok = false
            M_herm = zeros(ComplexF64, 2M, 2M)
            M_herm[1:M,     1:M]     .=  A
            M_herm[1:M,     M+1:2M]  .=  B_q
            M_herm[M+1:2M,  1:M]     .=  B_q'
            M_herm[M+1:2M,  M+1:2M]  .= -D_q
            M_herm .= 0.5 .* (M_herm .+ M_herm')
            for i in 1:2M; M_herm[i, i] += eta; end

            try
                K = cholesky(Hermitian(M_herm)).U
                sigma3 = Diagonal([fill(1.0, M); fill(-1.0, M)])
                W = K * sigma3 * K'
                W = Hermitian(0.5 .* (W .+ W'))
                F_W = eigen(W)

                pos_idx = findall(v -> v > 0, F_W.values)
                perm    = pos_idx[sortperm(F_W.values[pos_idx])]
                energies[iq]      = F_W.values[perm]
                wavefunctions[iq] = F_W.vectors[:, perm]
                cholesky_ok = true
                verbose && lock(log_lock) do
                    @printf("  q[%d/%d] = [%s] ... %s, M=%d, done (RPA/BdG-Cholesky)\n",
                        iq, Nq, join([@sprintf("%.4f", qi) for qi in q], ", "), d_method, M)
                end
            catch e
                if !isa(e, PosDefException)
                    rethrow(e)
                end
            end

            # ── Fallback: direct 2M×2M non-Hermitian diagonalization ──
            if !cholesky_ok
                RPA_mat = zeros(ComplexF64, 2M, 2M)
                RPA_mat[1:M,     1:M]     .=  A
                RPA_mat[1:M,     M+1:2M]  .=  B_q
                RPA_mat[M+1:2M,  1:M]     .= -B_q'
                RPA_mat[M+1:2M,  M+1:2M]  .=  D_q

                F_full = eigen(RPA_mat)
                all_evals = F_full.values
                max_imag  = maximum(abs.(imag.(all_evals)))
                pos_idx   = findall(v -> real(v) > 0, all_evals)
                perm      = pos_idx[sortperm(real.(all_evals[pos_idx]))]
                energies[iq]      = real.(all_evals[perm])
                wavefunctions[iq] = F_full.vectors[:, perm]
                verbose && lock(log_lock) do
                    @printf("  q[%d/%d] = [%s] ... %s, M=%d, done (RPA/direct 2M, max|imag|=%.2e, %d pos modes)\n",
                        iq, Nq, join([@sprintf("%.4f", qi) for qi in q], ", "), d_method,
                        M, max_imag, length(perm))
                end
            end
        end
    end

    return (
        qpoints       = qpoints,
        energies      = energies,
        wavefunctions = wavefunctions,
        triples       = triples_out,
        solver        = solver,
    )
end
