"""
Tests for excitations/particle_hole.jl — momentum lookup infrastructure

Covers _kindex_analytic: O(1) analytic BZ lookup with per-dimension grid sizes Ns.
Also tests the Ns_hf inference logic (fractional coordinate counting) and the
active_dims detection for embedded lower-dimensional grids.
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

# Build B_inv and Ns from avecs+box, mirroring solve_ph_excitations internals
function _setup_analytic(avecs, box)
    rvecs = _recip(avecs)
    kpts  = build_kpoints(avecs, box)
    B_inv = inv(hcat(rvecs...))
    Nk    = length(kpts)
    D     = length(kpts[1])
    frac_sets = [Set{Int}() for _ in 1:D]
    for k in kpts
        f = B_inv * k
        for d in 1:D
            push!(frac_sets[d], round(Int, f[d] * Nk))
        end
    end
    Ns = [length(s) for s in frac_sets]
    return kpts, B_inv, Ns, rvecs
end

function _test_analytic(avecs, box)
    kpts, B_inv, Ns, rvecs = _setup_analytic(avecs, box)
    Binv_full = inv(hcat(rvecs...))

    @testset "Ns inference: prod(Ns)==Nk" begin
        @test prod(Ns) == length(kpts)
    end

    @testset "exact grid match" begin
        for (i, k) in enumerate(kpts)
            @test MeanFieldTheories._kindex_analytic(k, B_inv, Ns) == i
        end
    end

    @testset "BZ translation invariance (±1 shift)" begin
        for (i, k) in enumerate(kpts), b in rvecs
            @test MeanFieldTheories._kindex_analytic(k .+ b, B_inv, Ns) == i
            @test MeanFieldTheories._kindex_analytic(k .- b, B_inv, Ns) == i
        end
    end

    @testset "k+q closure" begin
        for ki in eachindex(kpts), qi in eachindex(kpts)
            kpq = kpts[ki] .+ kpts[qi]
            idx = MeanFieldTheories._kindex_analytic(kpq, B_inv, Ns)
            diff = kpq .- kpts[idx]
            c    = Binv_full * diff
            @test norm(c .- round.(c)) < 1e-6
        end
    end
end

# ── Test suites ────────────────────────────────────────────────────────────────

@testset "_kindex_analytic" begin

    @testset "1D chain (a=1)" begin
        _test_analytic([[1.0]], (6,))
    end

    @testset "1D chain (a=2)" begin
        _test_analytic([[2.0]], (5,))
    end

    @testset "2D square (a=1)" begin
        _test_analytic([[1.0, 0.0], [0.0, 1.0]], (4, 4))
    end

    @testset "2D square, enlarged 2×2 magnetic cell" begin
        _test_analytic([[2.0, 0.0], [0.0, 2.0]], (3, 3))
    end

    @testset "2D rectangular (Nx≠Ny)" begin
        _test_analytic([[1.0, 0.0], [0.0, 1.0]], (6, 4))
    end

    @testset "2D triangular lattice" begin
        _test_analytic([[1.0, 0.0], [0.5, sqrt(3)/2]], (4, 4))
    end

    @testset "2D honeycomb Bravais" begin
        _test_analytic([[sqrt(3), 0.0], [sqrt(3)/2, 3/2]], (3, 3))
    end

    @testset "2D kagome Bravais (= triangular)" begin
        _test_analytic([[1.0, 0.0], [0.5, sqrt(3)/2]], (3, 3))
    end

    @testset "multi-BZ shift (n=2,3)" begin
        avecs = [[1.0, 0.0], [0.0, 1.0]]
        kpts, B_inv, Ns, rvecs = _setup_analytic(avecs, (4, 4))
        for (i, k) in enumerate(kpts), n in [2, 3, -2, -3]
            @test MeanFieldTheories._kindex_analytic(k .+ n .* rvecs[1], B_inv, Ns) == i
            @test MeanFieldTheories._kindex_analytic(k .+ n .* rvecs[2], B_inv, Ns) == i
        end
    end

    @testset "q outside first BZ → k+q closure" begin
        avecs = [[1.0, 0.0], [0.0, 1.0]]
        kpts, B_inv, Ns, rvecs = _setup_analytic(avecs, (4, 4))
        Binv_full = inv(hcat(rvecs...))
        for ki in eachindex(kpts), qi in eachindex(kpts)
            q_ext = kpts[qi] .+ rvecs[1]
            kpq   = kpts[ki] .+ q_ext
            idx   = MeanFieldTheories._kindex_analytic(kpq, B_inv, Ns)
            diff  = kpq .- kpts[idx]
            c     = Binv_full * diff
            @test norm(c .- round.(c)) < 1e-6
        end
    end

end

@testset "active_dims + Ns_hf inference (embedded 1D in 2D)" begin
    # 1D chain embedded in 2D: kpoints are [[kx, 0.0], ...], ky=0 always.
    # active_dims should detect only dim 1; Ns should be [N].
    avecs_2d = [[1.0, 0.0], [0.0, 1.0]]
    N = 8
    kpts = [[kx, 0.0] for kx in range(0, 2π, length=N+1)[1:N]]

    rvec_1d = [[2π, 0.0]]   # 1 reciprocal vector, 2D but only x-component active

    D_hf = length(kpts[1])   # = 2
    active_dims = [d for d in 1:D_hf if any(abs(k[d]) > 1e-12 for k in kpts)]
    n_active    = length(active_dims)

    @testset "detects 1 active dimension" begin
        @test n_active == 1
        @test active_dims == [1]
    end

    @testset "active_k extracts correct component" begin
        active_k(k) = length(k) == n_active ? collect(Float64, k) :
                      length(k) == D_hf     ? Float64[k[d] for d in active_dims] :
                      error("bad length")
        for k in kpts
            @test active_k(k) ≈ [k[1]]
            @test active_k([k[1]]) ≈ [k[1]]
        end
    end

    @testset "Ns_hf inference gives [N]" begin
        rv_active = [Float64[b[d] for d in active_dims] for b in rvec_1d]
        B_inv = inv(hcat(rv_active...))
        frac_sets = [Set{Int}() for _ in 1:n_active]
        for k in kpts
            f = B_inv * Float64[k[d] for d in active_dims]
            for d in 1:n_active
                push!(frac_sets[d], round(Int, f[d] * N))
            end
        end
        Ns = [length(s) for s in frac_sets]
        @test Ns == [N]
        @test prod(Ns) == N
    end

    @testset "_kindex_analytic with 1D active subspace" begin
        rv_active = [Float64[b[d] for d in active_dims] for b in rvec_1d]
        B_inv = inv(hcat(rv_active...))
        Ns    = [N]
        active_k(k) = Float64[k[min(length(k), d)] for d in active_dims]
        for (i, k) in enumerate(kpts)
            @test MeanFieldTheories._kindex_analytic(active_k(k), B_inv, Ns) == i
        end
        # BZ folding
        b_active = rv_active[1]
        for (i, k) in enumerate(kpts)
            ak = active_k(k)
            @test MeanFieldTheories._kindex_analytic(ak .+ b_active, B_inv, Ns) == i
            @test MeanFieldTheories._kindex_analytic(ak .- b_active, B_inv, Ns) == i
        end
    end
end
