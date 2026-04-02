"""
Square-lattice t-U-J-J' model at half-filling: magnon spectrum via RPA.

Hamiltonian:
  H = -t Σ_{⟨ij⟩,σ} c†_{iσ} c_{jσ}
    +  U Σ_i n_{i↑} n_{i↓}
    +  J  Σ_{⟨ij⟩}  S_i · S_j     (nearest-neighbor Heisenberg)
    + J' Σ_{⟨⟨ij⟩⟩} S_i · S_j     (next-nearest-neighbor Heisenberg)

where S_i · S_j = (1/2)(S⁺_i S⁻_j + S⁻_i S⁺_j) + Sᶻ_i Sᶻ_j.

The Heisenberg couplings J, J' > 0 favor antiferromagnetic order.
J' frustrates the Néel state and modifies the magnon dispersion.

Setup:
  - 2-sublattice (√2 × √2 rotated) magnetic unit cell.
  - Lattice vectors: a1 = [1,1], a2 = [1,-1].
  - 4 orbitals per cell: A↑, A↓, B↑, B↓.
  - Half-filling: 2 electrons per cell.

Run:
    julia --project=examples -t 8 examples/Heisenberg/run.jl
"""

using Printf
using LinearAlgebra
using MeanFieldTheories

# ── Parameters ────────────────────────────────────────────────────────────────
t   = 1.0       # hopping amplitude
U   = 8.0       # on-site Hubbard repulsion
J   = 1.0       # NN Heisenberg coupling (J > 0 = antiferromagnetic)
Jp  = 0.3 * J   # NNN Heisenberg coupling (J' > 0 = frustrating)

Nk_hf    = 12   # HF k-grid: Nk_hf × Nk_hf
n_q_ΓX   = 12   # q-path segment counts
n_q_XM   = 12
n_q_MΓ   = 12

# ── Magnetic unit cell (√2 × √2 rotated square lattice) ──────────────────────
a1 = [1.0,  1.0]
a2 = [1.0, -1.0]

unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [1.0, 0.0]];
    vectors = [a1, a2]
)

dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

nn_bonds     = bonds(unitcell, (:p, :p), 1)   # nearest-neighbor
nnn_bonds    = bonds(unitcell, (:p, :p), 2)   # next-nearest-neighbor
onsite_bonds = bonds(unitcell, (:p, :p), 0)   # on-site

# ── One-body: nearest-neighbor hopping ────────────────────────────────────────
onebody = generate_onebody(dofs, nn_bonds,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t : 0.0)

# ── Two-body interactions ─────────────────────────────────────────────────────

# On-site Hubbard U
hubbard = generate_twobody(dofs, onsite_bonds,
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin == qn2.spin) && (qn3.spin == qn4.spin) &&
        (qn1.spin !== qn3.spin) ? U/2 : 0.0,
    order = (cdag, :i, c, :i, cdag, :i, c, :i))

# S_i · S_j interaction generator for arbitrary bond set and coupling
function heisenberg_interaction(dofs, bond_set, coupling)
    generate_twobody(dofs, bond_set,
        (deltas, qn1, qn2, qn3, qn4) -> begin
            s1, s2, s3, s4 = qn1.spin, qn2.spin, qn3.spin, qn4.spin
            if (s1, s2, s3, s4) == (1, 2, 2, 1)         # S⁺_i S⁻_j
                return coupling / 2
            elseif (s1, s2, s3, s4) == (2, 1, 1, 2)     # S⁻_i S⁺_j
                return coupling / 2
            elseif s1 == s2 && s3 == s4                  # Sᶻ_i Sᶻ_j
                return s1 == s3 ? coupling / 4 : -coupling / 4
            else
                return 0.0
            end
        end,
        order = (cdag, :i, c, :i, cdag, :j, c, :j))
end

heisen_nn  = heisenberg_interaction(dofs, nn_bonds,  J)
heisen_nnn = heisenberg_interaction(dofs, nnn_bonds, Jp)

# Combine U + J + J'
twobody = (ops   = [hubbard.ops;    heisen_nn.ops;    heisen_nnn.ops],
           delta = [hubbard.delta;  heisen_nn.delta;  heisen_nnn.delta],
           irvec = [hubbard.irvec;  heisen_nn.irvec;  heisen_nnn.irvec])

# ── Reciprocal lattice vectors ────────────────────────────────────────────────
B_mat = 2π * inv(hcat(a1, a2))'
b1    = B_mat[:, 1]
b2    = B_mat[:, 2]
reciprocal_vecs = [b1, b2]

# ── k-grid & half-filling ─────────────────────────────────────────────────────
kpoints = build_kpoints([a1, a2], (Nk_hf, Nk_hf))
Nk      = length(kpoints)
n_elec  = 2 * Nk

# ── Hartree-Fock ──────────────────────────────────────────────────────────────
println("=" ^ 60)
println("Square-lattice t-U-J-J' model — magnon via RPA")
println(@sprintf("  t = %.2f,  U = %.2f,  J = %.2f,  J' = %.2f", t, U, J, Jp))
println(@sprintf("  Magnetic unit cell: a1=%s, a2=%s", a1, a2))
println(@sprintf("  HF k-grid: %d × %d = %d k-points, n_elec = %d", Nk_hf, Nk_hf, Nk, n_elec))
println("=" ^ 60)

hf = solve_hfk(dofs, onebody, twobody, kpoints, n_elec;
    n_restarts     = 10,
    field_strength = 1.0,
    n_warmup       = 10,
    tol            = 1e-12,
    verbose        = true)

# ── Check Neel order ──────────────────────────────────────────────────────────
mags = local_spin(dofs, hf.G_k)
print_spin(mags)

sA, sB = mags[1], mags[2]
m_stag = abs(sA.m + sB.m)

println(@sprintf("\nStaggered magnetization: %.6f", m_stag))
if m_stag < 0.1
    @error "Staggered magnetization too small — Neel order not established."
end

# ── q-path: Gamma -> X -> M -> Gamma (original BZ) ───────────────────────────
Gamma = [0.0, 0.0]
X_pt  = [π,   0.0]
M_pt  = [π,   π  ]

grid_path(p1, p2, n) = [p1 .+ (i / n) .* (p2 .- p1) for i in 0:n]

qpoints = [grid_path(Gamma, X_pt, n_q_ΓX);
           grid_path(X_pt,  M_pt, n_q_XM)[2:end];
           grid_path(M_pt, Gamma, n_q_MΓ)[2:end]]

d_ΓX = norm(X_pt .- Gamma)
d_XM = norm(M_pt .- X_pt)
d_MΓ = norm(Gamma .- M_pt)
arc = [collect(range(0.0,         d_ΓX;                length=n_q_ΓX+1));
       collect(range(d_ΓX,        d_ΓX + d_XM;         length=n_q_XM+1))[2:end];
       collect(range(d_ΓX + d_XM, d_ΓX + d_XM + d_MΓ; length=n_q_MΓ+1))[2:end]]

xtick_pos = [0.0, d_ΓX, d_ΓX + d_XM, d_ΓX + d_XM + d_MΓ]
xtick_lab = ["Γ", "X", "M", "Γ"]

# ── HF band structure ─────────────────────────────────────────────────────────
hf_bands, _ = energy_bands(dofs, onebody, twobody, kpoints, hf.G_k, qpoints)

open(joinpath(@__DIR__, "hf_bands.dat"), "w") do f
    norb_hf = size(hf_bands, 1)
    println(f, "# t-U-J-J' HF bands along Γ→X→M→Γ")
    println(f, "# t=$t  U=$U  J=$J  J'=$Jp")
    println(f, "# arc  E1 ... E$(norb_hf)")
    println(f, "# mu: $(hf.mu)")
    println(f, "# xtick_pos: $(xtick_pos)")
    println(f, "# xtick_lab: $(xtick_lab)")
    for qi in eachindex(qpoints)
        vals = join([@sprintf("%.8f", hf_bands[n, qi]) for n in 1:norb_hf], "  ")
        println(f, @sprintf("%.6f  %s", arc[qi], vals))
    end
end
println("Saved: hf_bands.dat")

# ── RPA magnon dispersion ─────────────────────────────────────────────────────
ph = solve_ph_excitations(dofs, onebody, twobody, hf, qpoints, reciprocal_vecs;
    solver = :RPA,
    verbose = true)

open(joinpath(@__DIR__, "magnon.dat"), "w") do f
    max_modes = maximum(length(ph.energies[qi]) for qi in eachindex(qpoints))
    println(f, "# t-U-J-J' magnon dispersion (RPA)")
    println(f, "# t=$t  U=$U  J=$J  J'=$Jp")
    println(f, "# arc  E1 ... E$(max_modes)  (sorted ascending)")
    println(f, "# xtick_pos: $(xtick_pos)")
    println(f, "# xtick_lab: $(xtick_lab)")
    for qi in eachindex(qpoints)
        all_E = sort(ph.energies[qi])
        vals  = join([@sprintf("%.8f", e) for e in all_E], "  ")
        println(f, @sprintf("%.6f  %s", arc[qi], vals))
    end
end
println("Saved: magnon.dat")

# ── Goldstone check ───────────────────────────────────────────────────────────
E_Gamma = sort(ph.energies[1])
println(@sprintf("\nLowest excitation at Γ: %.6e (should be ≈ 0 for Goldstone mode)", E_Gamma[1]))
println("Done.")
