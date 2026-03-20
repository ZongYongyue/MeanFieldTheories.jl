"""
Tests for excitations/particle_hole.jl

Covers _kpoint_index: BZ folding and grid lookup on several common lattices.
"""

using Test
using LinearAlgebra
using MeanFieldTheories

# Helper: reciprocal vectors from direct lattice vectors
# Returns [b1, b2, ...] satisfying a_i · b_j = 2π δ_{ij}
function _recip(unitcell_vecs)
    A = hcat(unitcell_vecs...)
    B = 2π * inv(A)'
    return [B[:, i] for i in 1:size(B, 2)]
end

# Core correctness check reused across all lattices:
#   1. every grid point maps to itself
#   2. k + b_i folds back to the same point (full BZ translation)
#   3. for every (ki, qj), k+q lands on a grid point and the result is consistent
#      with manual BZ folding
function _test_lattice(avecs, box)
    rvecs = _recip(avecs)
    B     = hcat(rvecs...)
    kpts  = build_kpoints(avecs, box)

    @testset "exact grid match" begin
        for (i, k) in enumerate(kpts)
            @test MeanFieldTheories._kpoint_index(kpts, k, rvecs) == i
        end
    end

    @testset "BZ translation invariance" begin
        for (i, k) in enumerate(kpts), b in rvecs
            @test MeanFieldTheories._kpoint_index(kpts, k .+ b, rvecs) == i
            @test MeanFieldTheories._kpoint_index(kpts, k .- b, rvecs) == i
        end
    end

    @testset "k+q closure: result differs from k+q by integer reciprocal shift" begin
        # kpts[idx] must equal kpq - Σ n_i b_i for some integers n_i.
        # Check this by verifying (kpq - kpts[idx]) expressed in reciprocal
        # coordinates is close to an integer vector.
        Binv = inv(B)
        for ki in eachindex(kpts), qi in eachindex(kpts)
            kpq = kpts[ki] .+ kpts[qi]
            idx = MeanFieldTheories._kpoint_index(kpts, kpq, rvecs)
            diff = kpq .- kpts[idx]          # should be integer combination of b_i
            c    = Binv * diff               # fractional reciprocal coordinates
            @test norm(c .- round.(c)) < 1e-6
        end
    end
end

@testset "_kpoint_index" begin

    @testset "1D chain" begin
        rvecs = _recip([[1.0]])
        kpts  = build_kpoints([[1.0]], (6,))
        @testset "exact grid match" begin
            for (i, k) in enumerate(kpts)
                @test MeanFieldTheories._kpoint_index(kpts, k, rvecs) == i
            end
        end
        @testset "BZ translation invariance" begin
            for (i, k) in enumerate(kpts)
                @test MeanFieldTheories._kpoint_index(kpts, k .+ rvecs[1], rvecs) == i
                @test MeanFieldTheories._kpoint_index(kpts, k .- rvecs[1], rvecs) == i
            end
        end
    end

    @testset "2D square lattice (a=1)" begin
        _test_lattice([[1.0, 0.0], [0.0, 1.0]], (4, 4))
    end

    @testset "2D square lattice, enlarged 2×2 magnetic cell" begin
        _test_lattice([[2.0, 0.0], [0.0, 2.0]], (3, 3))
    end

    @testset "2D triangular lattice" begin
        # a1=[1,0], a2=[1/2, √3/2]
        _test_lattice([[1.0, 0.0], [0.5, sqrt(3)/2]], (4, 4))
    end

    @testset "2D honeycomb lattice (primitive cell, 2 sites)" begin
        # Same Bravais lattice as triangular; primitive vectors:
        # a1=[√3, 0], a2=[√3/2, 3/2]  (nearest-neighbor distance = 1)
        _test_lattice([[sqrt(3), 0.0], [sqrt(3)/2, 3/2]], (3, 3))
    end

    @testset "2D kagome lattice (primitive cell = triangular)" begin
        # Kagome shares the triangular Bravais lattice with a 3-site basis.
        # The primitive direct vectors are the same as triangular:
        _test_lattice([[1.0, 0.0], [0.5, sqrt(3)/2]], (3, 3))
    end

    @testset "error on point not on grid" begin
        rvecs = _recip([[1.0, 0.0], [0.0, 1.0]])
        kpts  = build_kpoints([[1.0, 0.0], [0.0, 1.0]], (3, 3))
        @test_throws ErrorException MeanFieldTheories._kpoint_index(kpts, [0.1234, 0.5678], rvecs)
    end

    @testset "q outside first BZ (2nd and 3rd BZ)" begin
        # 2D square lattice, 4×4 grid; reciprocal vectors b1=[2π,0], b2=[0,2π]
        avecs = [[1.0, 0.0], [0.0, 1.0]]
        rvecs = _recip(avecs)
        kpts  = build_kpoints(avecs, (4, 4))
        # For every grid point k, shifting by n*b should map back to the same index
        for (i, k) in enumerate(kpts), n in [2, 3, -2, -3]
            @test MeanFieldTheories._kpoint_index(kpts, k .+ n .* rvecs[1], rvecs) == i
            @test MeanFieldTheories._kpoint_index(kpts, k .+ n .* rvecs[2], rvecs) == i
        end
        # k+q closure when q is in the 2nd BZ (q = k_j + b1, for arbitrary k_j)
        Binv = inv(hcat(rvecs...))
        for ki in eachindex(kpts), qi in eachindex(kpts)
            q_ext = kpts[qi] .+ rvecs[1]   # push q into 2nd BZ
            kpq   = kpts[ki] .+ q_ext
            idx   = MeanFieldTheories._kpoint_index(kpts, kpq, rvecs)
            diff  = kpq .- kpts[idx]
            c     = Binv * diff
            @test norm(c .- round.(c)) < 1e-6
        end
    end

end
