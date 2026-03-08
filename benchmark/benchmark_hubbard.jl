"""
Correctness benchmark for Hartree-Fock solvers.

  Test 1 — Hubbard model (t=1, U=8), 2D square lattice, 8-site √8×√8 cluster,
            half-filling: real-space HF (solve_hfr).

  Test 2 — V-model (t=1, V=4), 2D square lattice, 2×2 AFM magnetic unit cell,
            half-filling: momentum-space HF (solve_hfk).

Run from the repo root:
    julia --project benchmark/benchmark_hubbard.jl
"""

using Printf
using LinearAlgebra
using MeanFieldTheories

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: Hubbard model (t=1, U=8) — real-space HF, 8-site √8×√8 cluster
#
# H = -t Σ_{<ij>,σ} c†_{iσ}c_{jσ} + U Σ_i n_{i↑}n_{i↓},  t=1, U=8
# 8-site cluster with supercell vectors a0=(2,2), a1=(-2,2) (√8×√8, 45°-rotated)
# on the 2D square lattice.  Half-filling: 8 electrons, Sz=0.
# Reference total energy: -3.408
# ─────────────────────────────────────────────────────────────────────────────

println()
println("=" ^ 60)
println("Test 1: Hubbard (t=1, U=8) — real-space HF, 8-site √8×√8 cluster")
println("=" ^ 60)

const t_hub = 1.0
const U_hub = 8.0

dofs_r = SystemDofs([Dof(:site, 8), Dof(:spin, 2)], sortrule = [[2], 1])

# 8 sites at their positions in the √8×√8 fundamental domain.
# Supercell vectors [[2,2],[-2,2]] encode the tilted PBC; bonds() finds all NN pairs.
cluster_r = Lattice([Dof(:site, 8)],
                     [QN(site=i) for i in 1:8],
                     [[0.0,0.0], [-1.0,1.0], [0.0,1.0], [1.0,1.0],
                      [-1.0,2.0], [0.0,2.0], [1.0,2.0], [0.0,3.0]];
                     supercell_vectors=[[2.0,2.0],[-2.0,2.0]])

nn_bonds_r    = bonds(cluster_r, (:p, :p), 1)
onsite_bonds_r = bonds(cluster_r, (:p, :p), 0)

t_ops_r = generate_onebody(dofs_r, nn_bonds_r,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t_hub : 0.0,
    hc = true).ops

U_ops_r = generate_twobody(dofs_r, onsite_bonds_r,
    (deltas, qn1, qn2, qn3, qn4) ->
        (qn1.spin, qn2.spin, qn3.spin, qn4.spin) == (1, 1, 2, 2) ? U_hub : 0.0,
    order = (cdag, :i, c, :i, cdag, :i, c, :i)).ops

result_r = solve_hfr(dofs_r, vcat(t_ops_r, U_ops_r), [4, 4];
    temperature = 0.0, max_iter = 1000, n_restarts = 20,
    tol = 1e-8, mix_alpha = 0.5, verbose = true)

const E_ref_r = -3.408

println()
println(@sprintf("  real-space HF energy:  %+.10f", result_r.energies.total))
println(@sprintf("  Reference:             %+.10f", E_ref_r))
println(@sprintf("  |difference|:          %.3e", abs(result_r.energies.total - E_ref_r)))

@assert result_r.converged "Test 1 FAILED: solve_hfr did not converge"
@assert abs(result_r.energies.total - E_ref_r) < 1e-2 "Test 1 FAILED: energy mismatch"

println()
println("  ✓ Test 1 PASSED")

println()
println("=" ^ 60)
println("Test 2: V-model (t=1, V=4) — 2D square, 2×2 AFM cell, half-filling")
println("=" ^ 60)

# Model: H = -t Σ_{<ij>,σ} c†_{iσ}c_{jσ} + V Σ_{i≠j,σσ'} n_{iσ}n_{jσ'}
# t=1, V=4, 2D square lattice, half-filling (16 electrons on 4×4 sites).
# Magnetic unit cell: 2×2 (4 sites, 2 spins → d=8), allows Q=(π,π) AFM order.
# Reference energy per magnetic unit cell: -2.2626401 / 4

const t1 = 1.0
const V1 = 4.0

dofs_v = SystemDofs([Dof(:site, 4), Dof(:spin, 2, [:up, :dn])])

unitcell_v = Lattice([Dof(:site, 4)],
                     [QN(site=1), QN(site=2), QN(site=3), QN(site=4)],
                     [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
                     supercell_vectors=[[2.0, 0.0], [0.0, 2.0]])

# bonds() on the unit cell generates both intra-cell and inter-cell NN bonds.
# Total lattice (4×4 = 16 sites) is implicit in the 2×2 k-grid below.
nn_bonds_v = bonds(unitcell_v, (:p, :p), 1)

onebody_v = generate_onebody(dofs_v, nn_bonds_v,
    (delta, qn1, qn2) -> qn1.spin == qn2.spin ? -t1 : 0.0)

# CoulombInter: V Σ_{i≠j, σσ'} n_{iσ} n_{jσ'} (ordered-pair sum, both i→j and j→i).
# :i and :j are independent position indices; generate_twobody automatically enumerates
# all injective assignments of {:i,:j} to bond sites, covering both directions.
twobody_v = generate_twobody(dofs_v, nn_bonds_v,
    (deltas, qn1, qn2, qn3, qn4) ->
        qn1.spin == qn2.spin && qn3.spin == qn4.spin ? V1 : 0.0)

kpoints_v = build_kpoints([[2.0, 0.0], [0.0, 2.0]], (2, 2))
n_elec_v  = 16  # half-filling

result_v = solve_hfk(dofs_v, onebody_v, twobody_v, kpoints_v, n_elec_v;
    n_restarts=20, tol=1e-8, verbose=true)

# Reference: converged SCF energy per magnetic unit cell.
const E_ref_per_cell = -2.2626401 / 4

println()
println(@sprintf("  k-space HF energy (per cell): %+.10f", result_v.energies.total))
println(@sprintf("  Reference (per cell):         %+.10f", E_ref_per_cell))
println(@sprintf("  |difference|:                 %.3e",
                 abs(result_v.energies.total - E_ref_per_cell)))
println(@sprintf("  NCond (per cell):             %.8f  (expected %.8f)",
                 result_v.ncond, n_elec_v / 4))

@assert result_v.converged "Test 2 FAILED: solve_hfk did not converge"
@assert abs(result_v.ncond - n_elec_v / 4) < 1e-6 "Test 2 FAILED: wrong particle count"

println()
println("  ✓ Test 2 PASSED")

println()
println("=" ^ 60)
println("All benchmark tests passed ✓")
println("=" ^ 60)
