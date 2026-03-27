"""
Square-lattice Hubbard model at half-filling: Néel AFM magnon spectrum.

Validates solve_ph_excitations via the Goldstone theorem:
  - Néel order spontaneously breaks SU(2) spin symmetry.
  - The particle-hole excitation spectrum must contain a gapless
    Goldstone magnon at q = Γ = (0, 0) in the magnetic Brillouin zone.

Setup:
  - 2-sublattice (√2 × √2 rotated) magnetic unit cell.
  - Lattice vectors: a1 = [1,1], a2 = [1,-1].
  - 4 orbitals per cell: A↑, A↓, B↑, B↓.
  - Half-filling: 2 electrons per cell → bands 1,2 occupied, 3,4 empty.
  - Large U/t (default U = 8) → well-formed Néel state.

Goldstone check:
  - Run with n0_list = [1, 2] (all occupied bands as hole candidates).
  - n_list = [1, 2, 3, 4] (occ filter removes occupied states automatically).
  - The lowest eigenvalue of H_eff at q = Γ must be ≈ 0.

Run:
    julia --project=examples -t 8 examples/AFM_Magnon/run.jl
"""

using Printf
using LinearAlgebra
using MeanFieldTheories

# ── Parameters ────────────────────────────────────────────────────────────────
t        = 1.0
U        = 12.0     # large U: strong Néel order, magnon well-defined
Z        = -10
Nk_hf    = 12       # HF k-grid: Nk_hf × Nk_hf (small for debugging)

# q-path segment counts
n_q_ΓX   = 12
n_q_XM   = 12
n_q_MΓ   = 12

# ── Magnetic unit cell (√2 × √2 rotated square lattice) ──────────────────────
# Original square lattice with a=1; magnetic unit cell doubles the unit cell:
#   a1 = [1, 1],  a2 = [1, -1]
# Two sites per cell: A at (0,0), B at (1,0) in Cartesian.
a1 = [1.0,  1.0]
a2 = [1.0, -1.0]

unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [1.0, 0.0]];
    vectors = [a1, a2]
)

# System DOFs: 1 cell × 2 sublattices × 2 spins = 4 orbitals
dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

nn_bonds     = bonds(unitcell, (:p, :p), 1)
onsite_bonds = bonds(unitcell, (:p, :p), 0)

# ── One-body: spin-conserving nearest-neighbor hopping ────────────────────────
onebody = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0)

# ── Two-body: on-site Hubbard U n↑n↓ ─────────────────────────────────────────
twobody = generate_twobody(dofs, onsite_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin == qn2.spin) && (qn3.spin == qn4.spin) && (qn1.spin !== qn3.spin) ? U/2 : 0.0,
    order = (cdag, :i, c, :i, cdag, :i, c, :i))


# ── Reciprocal lattice vectors ────────────────────────────────────────────────
# b_i satisfying a_i · b_j = 2π δ_{ij}
B_mat = 2π * inv(hcat(a1, a2))'
b1    = B_mat[:, 1]   # = [π,  π]
b2    = B_mat[:, 2]   # = [π, -π]
reciprocal_vecs = [b1, b2]

# ── k-grids ───────────────────────────────────────────────────────────────────
kpoints    = build_kpoints([a1, a2], (Nk_hf,    Nk_hf))
Nk      = length(kpoints)
n_elec  = 2 * Nk    # half-filling: 2 electrons per 4-orbital cell

# ── Hartree-Fock ──────────────────────────────────────────────────────────────
println("=" ^ 60)
println("Square-lattice Hubbard model — Néel AFM magnon validation")
println(@sprintf("  t = %.1f,  U = %.1f,  U/t = %.1f", t, U, U/t))
println(@sprintf("  Magnetic unit cell: a1=%s, a2=%s", a1, a2))
println(@sprintf("  HF k-grid:    %d × %d = %d k-points,  n_elec = %d", Nk_hf, Nk_hf, Nk, n_elec))
println("=" ^ 60)

hf = solve_hfk(dofs, onebody, twobody, kpoints, n_elec;
    n_restarts     = 10,
    field_strength = 1.0,
    n_warmup       = 10,
    tol            = 1e-12,
    verbose        = true)

# ── Check Néel order ──────────────────────────────────────────────────────────
mags   = local_spin(dofs, hf.G_k)
print_spin(mags)

sA, sB = mags[1], mags[2]
m_stag = abs(sA.m + sB.m) 

println("\nLocal spin moments (A and B sublattices):")
println(@sprintf("  A: mx=%.4f  my=%.4f  mz=%.4f", sA.mx, sA.my, sA.mz))
println(@sprintf("  B: mx=%.4f  my=%.4f  mz=%.4f", sB.mx, sB.my, sB.mz))
println(@sprintf("  Staggered magnetization |mz_A - mz_B|/2 = %.6f", m_stag))

if m_stag < 0.1
    @error "Staggered magnetization is small — Néel order may not be established. Try a larger U."
end

# ── q-path for magnon dispersion ─────────────────────────────────────────────
# High-symmetry points of the ORIGINAL square lattice BZ (a=1):
#   Γ = (0,  0)
#   X = (π,  0)   [zone-face midpoint]
#   M = (π,  π)   [zone corner]
#
# In the magnetic unit cell, X is the M-point of the magnetic BZ, and
# M_orig = b1_mag folds back to magnetic Γ.  Points along X→M and near M
# lie in the 2nd magnetic BZ; _kpoint_index folds them back automatically.
#
# Divisibility w.r.t. HF grid (spacing b_i/Nk_hf):
#   Γ→X: step [π/n_q_ΓX, 0] → need Nk_hf % (2*n_q_ΓX) == 0
#   X→M: step [0, π/n_q_XM] → same condition
#   M→Γ: step [-π/n_q_MΓ, -π/n_q_MΓ] → need Nk_hf % n_q_MΓ == 0

Γ_pt = [0.0, 0.0]
X_pt = [π,   0.0]    # original BZ X
M_pt = [π,   π  ]    # original BZ M (= b1_mag, folds to magnetic Γ)

function grid_path(p1, p2, n)
    [p1 .+ (i / n) .* (p2 .- p1) for i in 0:n]
end

qpath_ΓX = grid_path(Γ_pt, X_pt, n_q_ΓX)
qpath_XM = grid_path(X_pt, M_pt, n_q_XM)
qpath_MΓ = grid_path(M_pt, Γ_pt, n_q_MΓ)
qpoints  = [qpath_ΓX; qpath_XM[2:end]; qpath_MΓ[2:end]]

# Arc-length parameterization for plotting
d_ΓX = norm(X_pt .- Γ_pt)
d_XM = norm(M_pt .- X_pt)
d_MΓ = norm(Γ_pt .- M_pt)
arc = [collect(range(0.0,           d_ΓX;                length=n_q_ΓX+1));
       collect(range(d_ΓX,          d_ΓX + d_XM;         length=n_q_XM+1))[2:end];
       collect(range(d_ΓX + d_XM,   d_ΓX + d_XM + d_MΓ; length=n_q_MΓ+1))[2:end]]

xtick_pos = [0.0, d_ΓX, d_ΓX + d_XM, d_ΓX + d_XM + d_MΓ]
xtick_lab = ["Γ", "X", "M", "Γ"]

# ── HF band structure along q-path ───────────────────────────────────────────
hf_bands, _ = energy_bands(dofs, onebody, twobody, kpoints, hf.G_k, qpoints)
# hf_bands: (norb, Nq), rows sorted ascending at each q

open(joinpath(@__DIR__, "hf_bands.dat"), "w") do f
    norb_hf = size(hf_bands, 1)
    println(f, "# HF mean-field bands along Γ→X→M→Γ (original BZ)")
    println(f, "# arc  E1 ... E$(norb_hf)")
    println(f, "# mu: $(hf.mu)")
    println(f, "# xtick_pos: $(xtick_pos)")
    println(f, "# xtick_lab: $(xtick_lab)")
    for qi in eachindex(qpoints)
        vals = join([@sprintf("%.8f", hf_bands[n, qi]) for n in 1:norb_hf], "  ")
        println(f, @sprintf("%.6f  %s", arc[qi], vals))
    end
end
println("Saved: examples/AFM_Magnon/hf_bands.dat")
println("Done.")

ph = solve_ph_excitations(dofs, onebody, twobody, hf, qpoints, reciprocal_vecs;
    solver = :RPA,
    verbose = true)


open(joinpath(@__DIR__, "magnon.dat"), "w") do f
    max_modes = maximum(length(ph.energies[qi]) for qi in eachindex(qpoints))
    println(f, "# arc  E1 ... E$(max_modes)  (sorted ascending)")
    println(f, "# xtick_pos: $(xtick_pos)")
    println(f, "# xtick_lab: $(xtick_lab)")
    for qi in eachindex(qpoints)
        all_E = sort(ph.energies[qi])
        vals  = join([@sprintf("%.8f", e) for e in all_E], "  ")
        println(f, @sprintf("%.6f  %s", arc[qi], vals))
    end
end
println("\nSaved: examples/AFM_Magnon/magnon.dat")


