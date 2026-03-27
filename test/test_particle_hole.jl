"""
Tests for excitations/particle_hole.jl — BZ folding and mesh lookup infrastructure

Covers _fold_to_BZ, _find_on_mesh, and _build_ph_pairs: the momentum-space
helpers that underpin the TDA/RPA excitation solver.
"""

using Test
using LinearAlgebra
using MeanFieldTheories

# ── Helpers ────────────────────────────────────────────────────────────────────

function _recip(avecs)
    A = hcat(avecs...)
    B = 2π * inv(A)'
    return [B[:, i] for i in 1:size(B, 2)]
end

"""Build B_mat, B_inv, kpoints, and kpoints_folded for a given lattice."""
function _setup_bz(avecs, box)
    rvecs = _recip(avecs)
    kpts  = build_kpoints(avecs, box)
    B_mat = hcat(rvecs...)
    B_inv = pinv(B_mat)
    kpts_folded = [MeanFieldTheories._fold_to_BZ(k, B_inv, B_mat) for k in kpts]
    return kpts, B_mat, B_inv, kpts_folded, rvecs
end

"""Test fold+mesh-hit for a given lattice: every k-point folds to itself,
BZ-shifted k folds back, and k+q always lands on a mesh point."""
function _test_fold_and_mesh(avecs, box)
    kpts, B_mat, B_inv, kpts_folded, rvecs = _setup_bz(avecs, box)

    @testset "self-consistency: fold(k) matches mesh" begin
        for (i, k) in enumerate(kpts)
            k_fold = MeanFieldTheories._fold_to_BZ(k, B_inv, B_mat)
            idx    = MeanFieldTheories._find_on_mesh(k_fold, kpts_folded)
            @test idx == i
        end
    end

    @testset "BZ translation invariance (±1 shift)" begin
        for (i, k) in enumerate(kpts), b in rvecs
            for sign in [1, -1]
                k_shift = k .+ sign .* b
                k_fold  = MeanFieldTheories._fold_to_BZ(k_shift, B_inv, B_mat)
                idx     = MeanFieldTheories._find_on_mesh(k_fold, kpts_folded)
                @test idx == i
            end
        end
    end

    @testset "k+q closure (all pairs land on mesh)" begin
        for ki in eachindex(kpts), qi in eachindex(kpts)
            kpq    = kpts[ki] .+ kpts[qi]
            k_fold = MeanFieldTheories._fold_to_BZ(kpq, B_inv, B_mat)
            idx    = MeanFieldTheories._find_on_mesh(k_fold, kpts_folded)
            @test idx !== nothing
        end
    end
end

# ── _fold_to_BZ and _find_on_mesh ─────────────────────────────────────────────

@testset "_fold_to_BZ + _find_on_mesh" begin

    @testset "1D chain (a=1)" begin
        _test_fold_and_mesh([[1.0]], (6,))
    end

    @testset "1D chain (a=2)" begin
        _test_fold_and_mesh([[2.0]], (5,))
    end

    @testset "2D square (a=1)" begin
        _test_fold_and_mesh([[1.0, 0.0], [0.0, 1.0]], (4, 4))
    end

    @testset "2D square, enlarged 2×2 magnetic cell" begin
        _test_fold_and_mesh([[2.0, 0.0], [0.0, 2.0]], (3, 3))
    end

    @testset "2D rectangular (Nx≠Ny)" begin
        _test_fold_and_mesh([[1.0, 0.0], [0.0, 1.0]], (6, 4))
    end

    @testset "2D triangular lattice" begin
        _test_fold_and_mesh([[1.0, 0.0], [0.5, sqrt(3)/2]], (4, 4))
    end

    @testset "2D honeycomb Bravais" begin
        _test_fold_and_mesh([[sqrt(3), 0.0], [sqrt(3)/2, 3/2]], (3, 3))
    end

    @testset "multi-BZ shift (n=2,3)" begin
        avecs = [[1.0, 0.0], [0.0, 1.0]]
        kpts, B_mat, B_inv, kpts_folded, rvecs = _setup_bz(avecs, (4, 4))
        for (i, k) in enumerate(kpts), n in [2, 3, -2, -3]
            for b in rvecs
                k_shift = k .+ n .* b
                k_fold  = MeanFieldTheories._fold_to_BZ(k_shift, B_inv, B_mat)
                idx     = MeanFieldTheories._find_on_mesh(k_fold, kpts_folded)
                @test idx == i
            end
        end
    end

    @testset "pinv handles 1D embedded in 2D" begin
        # 1D chain along x embedded in 2D: reciprocal vector is (2π, 0)
        rvecs = [[2π, 0.0]]
        B_mat = hcat(rvecs...)          # 2×1 matrix
        B_inv = pinv(B_mat)             # 1×2 matrix
        N = 6
        kpts = [[2π * i / N, 0.0] for i in 0:N-1]
        kpts_folded = [MeanFieldTheories._fold_to_BZ(k, B_inv, B_mat) for k in kpts]

        for (i, k) in enumerate(kpts)
            # Shift by reciprocal vector
            k_shift = k .+ rvecs[1]
            k_fold  = MeanFieldTheories._fold_to_BZ(k_shift, B_inv, B_mat)
            idx     = MeanFieldTheories._find_on_mesh(k_fold, kpts_folded)
            @test idx == i
        end
    end
end

# ── _build_ph_pairs ───────────────────────────────────────────────────────────

@testset "_build_ph_pairs" begin

    @testset "half-filled 2-band model" begin
        # 2 bands, 4 k-points; band 1 occupied, band 2 unoccupied
        norb = 2; Nk = 4
        evals_k  = zeros(norb, Nk)
        evals_kq = zeros(norb, Nk)
        for ki in 1:Nk
            evals_k[1, ki]  = -1.0    # occupied
            evals_k[2, ki]  =  1.0    # unoccupied
            evals_kq[1, ki] = -0.5    # occupied at k+q
            evals_kq[2, ki] =  1.5    # unoccupied at k+q
        end
        mu = 0.0

        triples = MeanFieldTheories._build_ph_pairs(evals_k, evals_kq, mu)
        # Each k gives 1 hole (band 1) × 1 particle (band 2) = 1 pair
        @test length(triples) == Nk
        for (ki, n0, n) in triples
            @test n0 == 1    # hole band
            @test n  == 2    # particle band
        end
    end

    @testset "fully occupied → no pairs" begin
        norb = 2; Nk = 3
        evals_k  = fill(-1.0, norb, Nk)   # all occupied
        evals_kq = fill(-0.5, norb, Nk)   # all occupied at k+q
        mu = 0.0

        triples = MeanFieldTheories._build_ph_pairs(evals_k, evals_kq, mu)
        @test isempty(triples)
    end

    @testset "n_list restricts particle bands" begin
        norb = 3; Nk = 2
        evals_k  = zeros(norb, Nk)
        evals_kq = zeros(norb, Nk)
        for ki in 1:Nk
            evals_k[1, ki]  = -1.0   # occupied
            evals_k[2, ki]  =  1.0   # unoccupied
            evals_k[3, ki]  =  2.0   # unoccupied
            evals_kq[1, ki] = -0.5   # occupied at k+q
            evals_kq[2, ki] =  1.5   # unoccupied at k+q
            evals_kq[3, ki] =  3.0   # unoccupied at k+q
        end
        mu = 0.0

        # Without n_list: 1 hole × 2 particles = 2 pairs per k
        triples_all = MeanFieldTheories._build_ph_pairs(evals_k, evals_kq, mu)
        @test length(triples_all) == 2 * Nk

        # With n_list=[2]: only band 2 allowed as particle
        triples_restricted = MeanFieldTheories._build_ph_pairs(evals_k, evals_kq, mu;
                                                                n_list=[2])
        @test length(triples_restricted) == Nk
        for (ki, n0, n) in triples_restricted
            @test n == 2
        end
    end

    @testset "k-dependent occupation" begin
        # Band 1 crosses Fermi level: occupied at k=1, unoccupied at k=2
        norb = 2; Nk = 2
        evals_k  = [-0.5  0.5;     # band 1: occ at k1, unocc at k2
                      1.0  1.0]     # band 2: always unocc
        evals_kq = [ 0.5  0.5;     # band 1 at k+q: always unocc
                      1.5  1.5]     # band 2 at k+q: always unocc
        mu = 0.0

        triples = MeanFieldTheories._build_ph_pairs(evals_k, evals_kq, mu)
        # k1: hole=band1, particle=band1,2 → 2 pairs
        # k2: no holes → 0 pairs
        @test length(triples) == 2
        for (ki, n0, n) in triples
            @test ki == 1
            @test n0 == 1
        end
    end
end
