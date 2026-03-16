"""
Ground state analysis utilities: local observables from HF density matrices.
"""

# ──────────────── Local magnetic moment ────────────────

"""
    local_spin(dofs, G; spin_dof=:spin, kwargs...) -> Vector{NamedTuple}

Compute the local magnetic moment on each site from a converged HF density matrix,
with optional label-based filtering via keyword arguments.

A **site** is defined by all quantum numbers *except* the spin DOF, so the function
works for any system (simple site+spin, multi-orbital, supercell, …).

Compatible with both HF solvers:
- **HFr**: pass `result.G`   — a real-space density matrix of shape `(N, N)`
- **HFk**: pass `result.G_k` — a k-space density matrix of shape `(d, d, Nk)`;
  the BZ average `G_loc = Σ_k G_k / Nk` is computed automatically.

# Requirements
`dofs` must contain exactly one spin DOF of size 2 (default name `:spin`).
Spin label convention: value `1` → spin-up, value `2` → spin-down.

# Keyword arguments
- `spin_dof`: name of the spin DOF (default `:spin`)
- any other keyword: filter criterion on the site label. The value can be a
  single integer or a vector of integers. Only sites matching **all** criteria
  are returned.

# Returns
A `Vector` of `NamedTuple`s, one per (matching) site, ordered by each site's
first linear index in `dofs.valid_states`. Each entry contains:

| field       | description                                              |
|-------------|----------------------------------------------------------|
| `label`     | `NamedTuple` of non-spin QNs, e.g. `(site=2,)` or `(cell=1, sub=2)` |
| `n`         | total occupation `⟨n↑⟩ + ⟨n↓⟩`                         |
| `mx,my,mz`  | spin vector `(⟨Sˣ⟩, ⟨Sʸ⟩, ⟨Sᶻ⟩)`                       |
| `m`         | magnitude `√(mx²+my²+mz²)`                              |
| `theta_deg` | polar angle from +z axis (degrees, 0°=up, 180°=down)    |
| `phi_deg`   | azimuthal angle in xy-plane from +x (degrees)           |

# Spin-component formulas
Given `G_loc[a,b] = ⟨c†_a c_b⟩`:
- `mz = (G[↑,↑] − G[↓,↓]) / 2`
- `mx = Re(G[↑,↓])`
- `my = Im(G[↑,↓])`

# Examples
```julia
result = solve_hfk(...)

# All sites
mags = local_spin(dofs, result.G_k)

# Only orbital=1 (a orbital) across all sites
mags_a = local_spin(dofs, result.G_k; orbital=1)

# Only B sublattice
mags_B = local_spin(dofs, result.G_k; sub=2)

# B sublattice at cells 3 and 4
mags_B34 = local_spin(dofs, result.G_k; sub=2, cell=[3,4])

print_spin(mags_B)
```
"""
function local_spin(dofs::SystemDofs, G::AbstractArray;
                            spin_dof::Symbol = :spin, kwargs...)

    # ── 1. Locate and validate spin DOF ────────────────────────────────────
    spin_pos = findfirst(d -> d.name == spin_dof, dofs.dofs)
    isnothing(spin_pos) &&
        error("No DOF named :$spin_dof in dofs. Available: $(getfield.(dofs.dofs, :name))")
    dofs.dofs[spin_pos].size == 2 ||
        error("Spin DOF :$spin_dof must have exactly 2 states (up/down), got $(dofs.dofs[spin_pos].size)")

    # ── 2. Build local density matrix ──────────────────────────────────────
    # HFk: G is (d, d, Nk) → BZ-average to (d, d)
    # HFr: G is (N, N)     → use directly
    G_loc = if ndims(G) == 3
        dropdims(sum(G; dims=3); dims=3) ./ size(G, 3)
    elseif ndims(G) == 2
        G
    else
        error("G must be a 2D (HFr) or 3D (HFk) array, got $(ndims(G))D")
    end

    # ── 3. Group (up, down) index pairs by site label ───────────────────────
    site_names = Tuple(d.name for d in dofs.dofs if d.name != spin_dof)

    # site_label (NamedTuple) → Dict(spin_value::Int → linear_index::Int)
    site_map = Dict{Any, Dict{Int,Int}}()
    for qn in dofs.valid_states
        site_key = NamedTuple{site_names}(Tuple(qn[n] for n in site_names))
        spin_val = qn[spin_dof]
        get!(site_map, site_key, Dict{Int,Int}())[spin_val] = dofs.qn_to_idx[qn]
    end

    # ── 4. Build label filter from remaining kwargs ─────────────────────────
    label_filter = pairs(kwargs)   # e.g. (sub=2, cell=[3,4])
    for key in keys(kwargs)
        key in site_names ||
            error("Filter key :$key is not a site label DOF. Available: $(site_names)")
    end

    # ── 5. Compute spin moments, sorted by first linear index ───────────────
    sorted_sites = sort(collect(site_map); by = kv -> minimum(values(kv[2])))

    results = NamedTuple[]
    for (site_key, spin_idx) in sorted_sites
        length(spin_idx) == 2 || continue   # skip if up or down index missing

        # Apply label filter (skip if any criterion is not satisfied)
        all(label_filter) do (key, val)
            haskey(site_key, key) || return false
            v = site_key[key]
            val isa AbstractVector ? v in val : v == val
        end || continue

        i_up = spin_idx[1]   # spin=1 → up
        i_dn = spin_idx[2]   # spin=2 → down

        n_up = real(G_loc[i_up, i_up])
        n_dn = real(G_loc[i_dn, i_dn])
        G_ud = G_loc[i_up, i_dn]   # ⟨c†_↑ c_↓⟩

        n  = n_up + n_dn
        mz = (n_up - n_dn) / 2
        mx = real(G_ud)
        my = imag(G_ud)
        m  = sqrt(mx^2 + my^2 + mz^2)

        # Polar angles (degenerate at m≈0 or on z-axis → set to 0)
        theta = m > 1e-10 ? atan(sqrt(mx^2 + my^2), mz) : 0.0
        phi   = (mx^2 + my^2) > 1e-20 ? atan(my, mx) : 0.0

        push!(results, (
            label     = site_key,
            n         = n,
            mx        = mx,
            my        = my,
            mz        = mz,
            m         = m,
            theta_deg = rad2deg(theta),
            phi_deg   = rad2deg(phi),
        ))
    end

    return results
end

# ──────────────── Pretty printer ────────────────

"""
    print_spin(mags; digits=4, io=stdout)

Pretty-print the output of [`local_spin`](@ref).

Columns: site label | n | mx | my | mz | |m| | θ(°) | φ(°)

# Example output
```
site                     n       mx       my       mz      |m|     θ(°)     φ(°)
────────────────────────────────────────────────────────────────────────────────
1                   1.0000   0.0000   0.0000   0.2341   0.2341     0.00     0.00
2                   1.0000   0.0000   0.0000  -0.2341   0.2341   180.00     0.00
```
"""
function print_spin(mags::Vector; digits::Int = 4, io::IO = stdout)
    isempty(mags) && (println(io, "(no sites)"); return)

    label_keys = keys(first(mags).label)
    site_header = join(string.(label_keys), "/")

    col_w = max(20, length(site_header) + 2)
    header = @sprintf("%-*s  %8s  %8s  %8s  %8s  %8s  %8s  %8s",
                      col_w, site_header, "n", "mx", "my", "mz", "|m|", "θ(°)", "φ(°)")
    println(io, header)
    println(io, "─" ^ length(header))

    for mag in mags
        site_str = join([string(mag.label[k]) for k in label_keys], "/")
        println(io, @sprintf("%-*s  %8.*f  %8.*f  %8.*f  %8.*f  %8.*f  %8.*f  %8.*f",
                              col_w, site_str,
                              digits, mag.n,
                              digits, mag.mx,
                              digits, mag.my,
                              digits, mag.mz,
                              digits, mag.m,
                              2,      mag.theta_deg,
                              2,      mag.phi_deg))
    end
end

# ──────────────── Magnetic structure factor & ordering wavevector ────────────────

"""
    spin_structure_factor(mags, lattice, qpoints) -> Vector{Float64}

Compute the magnetic structure factor

```math
S(\\mathbf{q}) = \\frac{1}{N_s} \\left|\\sum_s e^{i\\mathbf{q}\\cdot\\mathbf{d}_s}\\, \\vec{m}_s\\right|^2
```

at each q-point, where the sum runs over all sites, `\\mathbf{d}_s` is the
real-space coordinate of site `s` (looked up from `lattice`), and
`\\vec{m}_s = (m_x, m_y, m_z)` is its local magnetic moment
(e.g. from [`local_spin`](@ref)).

# Arguments
- `mags`: output of `local_spin`
- `lattice`: the `Lattice` used to build the system; site coordinates are
  extracted automatically from `mag.label` via `get_coordinate`
- `qpoints`: collection of q-vectors to evaluate at (same format as `build_kpoints`)

# Returns
`Vector{Float64}` of length `length(qpoints)`.

# Example
```julia
mags = local_spin(dofs, result.G_k)
qpts = build_kpoints([a1, a2], (100, 100))
Sq   = spin_structure_factor(mags, lattice, qpts)
```
"""
function spin_structure_factor(mags, lattice::Lattice, qpoints)
    coords = [get_coordinate(lattice, QuantumNumber(mag.label)) for mag in mags]
    Ns     = length(mags)
    Sq     = zeros(Float64, length(qpoints))
    for (iq, q) in enumerate(qpoints)
        Fx = Fy = Fz = zero(ComplexF64)
        for (mag, pos) in zip(mags, coords)
            phase = exp(im * sum(q .* pos))
            Fx += phase * mag.mx
            Fy += phase * mag.my
            Fz += phase * mag.mz
        end
        Sq[iq] = (abs2(Fx) + abs2(Fy) + abs2(Fz)) / Ns
    end
    return Sq
end

"""
    magnetic_ordering_wavevector(mags, lattice, qpoints) -> NamedTuple

Find the magnetic ordering wavevector **Q** as the q-point that maximises
[`spin_structure_factor`](@ref).

# Arguments
- `mags`: output of `local_spin`
- `lattice`: the `Lattice` used to build the system
- `qpoints`: collection of q-vectors to search over (same format as `build_kpoints`)

# Returns
A `NamedTuple` with fields:

| field     | description                                         |
|-----------|-----------------------------------------------------|
| `Q`       | the q-vector with the largest `S(q)`               |
| `Sq_max`  | the peak value `S(Q)`                              |
| `Sq`      | full `Vector{Float64}` of `S(q)` at all q-points  |

# Example
```julia
mags = local_spin(dofs, result.G_k)
qpts = build_kpoints([a1, a2], (100, 100))
res  = magnetic_ordering_wavevector(mags, lattice, qpts)
println("Q = ", res.Q, "  S(Q) = ", res.Sq_max)
```
"""
function magnetic_ordering_wavevector(mags, lattice::Lattice, qpoints)
    Sq    = spin_structure_factor(mags, lattice, qpoints)
    i_max = argmax(Sq)
    return (Q = qpoints[i_max], Sq_max = Sq[i_max], Sq = Sq)
end

# ──────────────── Local charge ────────────────

"""
    local_charge(dofs, G; kwargs...) -> Vector{NamedTuple}

Compute the local charge (occupation) for each valid quantum state from a
converged HF density matrix.

All DOFs (spin, orbital, sublattice, …) are treated equally: no DOF is singled
out or summed over.  Each valid state in `dofs` gets its own occupation
`n = G_loc[i, i]`.  The full quantum-number label is preserved so that
[`charge_structure_factor`](@ref) can extract spatial coordinates via the
`Lattice`.

Compatible with both HF solvers:
- **HFk**: pass `result.G_k` — BZ average is computed automatically.
- **HFr**: pass `result.G` — used directly.

# Keyword arguments
Any keyword filters on the quantum-number label (same convention as
[`local_spin`](@ref)).

# Returns
`Vector` of `NamedTuple`s with fields `label` and `n`, one per matching state,
sorted by linear index.

# Example
```julia
charges  = local_charge(dofs, result.G_k)
charges_A = local_charge(dofs, result.G_k; sub=1)
```
"""
function local_charge(dofs::SystemDofs, G::AbstractArray; kwargs...)

    G_loc = if ndims(G) == 3
        dropdims(sum(G; dims=3); dims=3) ./ size(G, 3)
    elseif ndims(G) == 2
        G
    else
        error("G must be a 2D (HFr) or 3D (HFk) array, got $(ndims(G))D")
    end

    all_names    = Tuple(d.name for d in dofs.dofs)
    label_filter = pairs(kwargs)
    for key in keys(kwargs)
        key in all_names ||
            error("Filter key :$key is not a DOF. Available: $(all_names)")
    end

    results = NamedTuple[]
    for qn in sort(dofs.valid_states; by = q -> dofs.qn_to_idx[q])
        label = NamedTuple{all_names}(Tuple(qn[n] for n in all_names))
        all(label_filter) do (key, val)
            v = label[key]
            val isa AbstractVector ? v in val : v == val
        end || continue
        idx = dofs.qn_to_idx[qn]
        push!(results, (label = label, n = real(G_loc[idx, idx])))
    end
    return results
end

# ──────────────── Charge structure factor & ordering wavevector ────────────────

"""
    charge_structure_factor(charges, lattice, qpoints) -> Vector{Float64}

Compute the charge structure factor

```math
N(\\mathbf{q}) = \\frac{1}{N_s} \\left|\\sum_s e^{i\\mathbf{q}\\cdot\\mathbf{d}_s}\\, \\delta n_s\\right|^2
```

where ``\\delta n_s = n_s - \\bar{n}`` is the deviation of the local charge from
its spatial average.  Subtracting the mean suppresses the trivial ``\\mathbf{q}=0``
peak from uniform filling, revealing true CDW modulations.

# Arguments
- `charges`: output of [`local_charge`](@ref)
- `lattice`: the `Lattice` used to build the system; site coordinates are
  extracted automatically via `get_coordinate`
- `qpoints`: collection of q-vectors to evaluate at (same format as `build_kpoints`)

# Returns
`Vector{Float64}` of length `length(qpoints)`.

# Example
```julia
charges = local_charge(dofs, result.G_k)
qpts    = build_kpoints([a1, a2], (100, 100))
Nq      = charge_structure_factor(charges, lattice, qpts)
```
"""
function charge_structure_factor(charges, lattice::Lattice, qpoints)
    coords = [get_coordinate(lattice, QuantumNumber(c.label)) for c in charges]
    n_vals = [c.n for c in charges]
    n_mean = sum(n_vals) / length(n_vals)
    δn     = n_vals .- n_mean

    Ns = length(charges)
    Nq = zeros(Float64, length(qpoints))
    for (iq, q) in enumerate(qpoints)
        F = zero(ComplexF64)
        for (dn, pos) in zip(δn, coords)
            F += exp(im * sum(q .* pos)) * dn
        end
        Nq[iq] = abs2(F) / Ns
    end
    return Nq
end

"""
    charge_ordering_wavevector(charges, lattice, qpoints) -> NamedTuple

Find the charge ordering wavevector **Q** as the q-point that maximises
[`charge_structure_factor`](@ref).

# Returns
A `NamedTuple` with fields:

| field     | description                                          |
|-----------|------------------------------------------------------|
| `Q`       | the q-vector with the largest `N(q)`                |
| `Nq_max`  | the peak value `N(Q)`                               |
| `Nq`      | full `Vector{Float64}` of `N(q)` at all q-points   |

# Example
```julia
charges = local_charge(dofs, result.G_k)
qpts    = build_kpoints([a1, a2], (100, 100))
res     = charge_ordering_wavevector(charges, lattice, qpts)
println("Q = ", res.Q, "  N(Q) = ", res.Nq_max)
```
"""
function charge_ordering_wavevector(charges, lattice::Lattice, qpoints)
    Nq    = charge_structure_factor(charges, lattice, qpoints)
    i_max = argmax(Nq)
    return (Q = qpoints[i_max], Nq_max = Nq[i_max], Nq = Nq)
end

# ──────────────── Plotting stub (implemented in ext/MakieExt.jl) ────────────────

"""
    plot_spin(mags, positions; kwargs...) -> Figure

Visualize local magnetic moments as arrows on a lattice.

Requires Makie to be loaded before calling:
```julia
using GLMakie    # interactive 3D
using CairoMakie # save to file
using WGLMakie   # Jupyter notebook
```

# Arguments
- `mags`: output of [`local_spin`](@ref)
- `positions`: coordinates of each site, same order as `mags`.
  Each element is a 2- or 3-component vector `[x, y]` or `[x, y, z]`.

# Keyword arguments
| kwarg             | default     | description                                                    |
|-------------------|-------------|----------------------------------------------------------------|
| `arrow_color`     | `:crimson`  | arrow color                                                    |
| `arrow_frac`      | `0.33`      | max arrow length = `min_dist × arrow_frac`                     |
| `shaft_lw`        | `2.0`       | arrow shaft line width                                         |
| `head_px`         | `13`        | arrowhead marker size (pixels)                                 |
| `head_frac`       | `0.32`      | fraction of arrow length taken by the head                     |
| `site_markersize` | `20`        | outer circle marker size                                       |
| `dot_size_frac`   | `0.45`      | inner ⊙/⊗ max size as fraction of `site_markersize`           |
| `bonds`           | `nothing`   | bond list from [`bonds`](@ref) — passed directly, no helper needed |
| `bond_color`      | `:gray60`   | bond line color                                                |
| `bond_width`      | `1.5`       | bond line width                                                |
| `unitcell_vecs`   | `nothing`   | two lattice vectors to draw unit cell outline                  |
| `unitcell_origin` | `nothing`   | origin of unit cell (defaults to first site)                   |
| `axis_padding`    | `0.5`       | extra space around sites on each axis                          |
| `title`           | `""`        | figure title                                                   |

# Layout
Single top-view (xy) panel. All three spin components are encoded simultaneously:
- **Arrow**: direction = `(mx, my)`, length ∝ `sqrt(mx²+my²)` (in-plane projection)
- **⊙ dot**: `mz > 0` (spin out of page), size ∝ `mz`
- **⊗ cross**: `mz < 0` (spin into page), size ∝ `|mz|`

Sites with purely in-plane spins show only arrows; purely z-polarized sites show only
⊙/⊗; canted spins show both simultaneously.

# Example
```julia
using CairoMakie
result = solve_hfk(...)
mags   = local_spin(dofs, result.G_k)

nn = bonds(lattice, (:p, :p), 1)
fig = plot_spin(mags, unitcell.coordinates;
    bonds         = nn,
    unitcell_vecs = [[3/2, √3/2], [0.0, √3]],
    title         = "KMH ground state")
save("magnetization.pdf", fig)
```
"""
function plot_spin end
