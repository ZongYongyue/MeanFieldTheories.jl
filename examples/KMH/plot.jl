"""
Plot KMH results (Rachel & Le Hur):
  - Fig. 5: KM non-interacting band structure (λ = 0.05, 0.2, 0.5, 1.0)
  - Fig. 6: SDW phase boundary Uc vs λ  (TODO: uncomment when res.dat is ready)

Reads:
  examples/KMH/bands.dat     — non-interacting bands
  examples/KMH/res.dat       — λ ∈ {0.2..1.0}, columns: lambda U m_neel_xy converged
  examples/SM_AFM/res.dat    — λ = 0 (pure Hubbard), columns: U m_AF E_gs phase

Run:
    julia --project=examples examples/KMH/plot.jl
"""

using CairoMakie
using LinearAlgebra
using Printf

# ── Tick positions (must match run.jl arc computation) ───────────────────────
# Reciprocal lattice for a1=[0,√3], a2=[3/2,√3/2]
const a1 = [0.0, √3]; const a2 = [3/2, √3/2]
const B_mat = 2π * inv(hcat(a1, a2))'
const b1 = B_mat[:, 1]; const b2 = B_mat[:, 2]
const Γ  = [0.0, 0.0]
const K  = (2b1 + b2) / 3
const Kp = (b1 + 2b2) / 3
const M  = (b1 + b2) / 2

const xtick_pos = cumsum([0.0, norm(Γ-K), norm(Kp-Γ), norm(K-Kp), norm(M-K), norm(Γ-M)])
const xtick_lab = ["K", "Γ", "K'", "K", "M", "Γ"]

# ── Read bands.dat ────────────────────────────────────────────────────────────
bands_data = Dict{Float64, Tuple{Vector{Float64}, Matrix{Float64}}}()

open(joinpath(@__DIR__, "bands.dat")) do f
    for line in eachline(f)
        startswith(line, "#") && continue
        isempty(strip(line)) && continue
        p   = split(line)
        λ   = parse(Float64, p[1])
        arc = parse(Float64, p[2])
        Es  = [parse(Float64, p[i]) for i in 3:lastindex(p)]
        if !haskey(bands_data, λ)
            bands_data[λ] = (Float64[], zeros(length(Es), 0))
        end
        arcs, bmat = bands_data[λ]
        push!(arcs, arc)
        bands_data[λ] = (arcs, hcat(bmat, Es))  # Es is a Vector, hcat appends as column
    end
end

λ_list = sort(collect(keys(bands_data)))

# ── Plot Fig. 5: band structure ───────────────────────────────────────────────
fig = Figure(size=(600, 900))
for (row, λ) in enumerate(λ_list)
    arcs, bmat = bands_data[λ]
    nbands = size(bmat, 1)
    half   = nbands ÷ 2

    ax = Axis(fig[row, 1];
        ylabel       = "E / t",
        title        = "λ = $λ",
        xticks       = (xtick_pos, xtick_lab),
        xgridvisible = false,
        ygridvisible = false,
        limits       = ((xtick_pos[1], xtick_pos[end]), nothing))
    row < length(λ_list) && hidexdecorations!(ax; ticks=false, grid=false)

    hlines!(ax, [0.0]; color=:gray, linestyle=:dash, linewidth=0.8)
    for x in xtick_pos[2:end-1]
        vlines!(ax, [x]; color=:gray, linestyle=:dot, linewidth=0.8)
    end
    for n in 1:half
        lines!(ax, arcs, bmat[n, :];      color=:steelblue, linewidth=1.5, overdraw=true)
    end
    for n in half+1:nbands
        lines!(ax, arcs, bmat[n, :];      color=:crimson,   linewidth=1.5, overdraw=true)
    end
end

display(fig)
save(joinpath(@__DIR__, "../../docs/src/fig/km_bands.png"), fig)
println("Saved: docs/src/fig/km_bands.png")

# ── Plot Fig. 6: phase boundary ───────────────────────────────────────────────
const mxy_threshold = 0.01

kmh_data = Dict{Float64, Vector{@NamedTuple{U::Float64, mxy::Float64}}}()
open(joinpath(@__DIR__, "res.dat")) do f
    for line in eachline(f)
        startswith(line, "#") && continue
        isempty(strip(line)) && continue
        p = split(line)
        λ, U, mxy = parse(Float64, p[1]), parse(Float64, p[2]), parse(Float64, p[3])
        push!(get!(kmh_data, λ, []), (U=U, mxy=mxy))
    end
end

# ── Plot order parameter vs U ─────────────────────────────────────────────────
λ_vals = sort(collect(keys(kmh_data)))

markers_op = [:circle, :rect, :utriangle, :diamond, :pentagon, :star5]
colors_op  = [:black, :steelblue, :forestgreen, :darkorange, :crimson, :mediumpurple]

fig_op = Figure(size=(560, 420))
ax_op  = Axis(fig_op[1, 1];
    xlabel       = "U / t",
    ylabel       = "m_Néel (xy)",
    title        = "KMH SDW order parameter",
    xgridvisible = false,
    ygridvisible = false,
    limits       = (nothing, (-0.02, nothing)))

for (i, λ) in enumerate(λ_vals)
    rows = sort(kmh_data[λ], by = r -> r.U)
    Us   = [r.U   for r in rows]
    mxys = [r.mxy for r in rows]
    scatterlines!(ax_op, Us, mxys;
        marker    = markers_op[i],
        color     = (colors_op[i], 0.85),
        linewidth = 1.2,
        markersize = 10,
        label     = "λ = $λ",
        overdraw  = true)
end
axislegend(ax_op; position=:lt)

display(fig_op)
save(joinpath(@__DIR__, "../../docs/src/fig/kmh_order_parameter.png"), fig_op)
println("Saved: docs/src/fig/kmh_order_parameter.png")

# ── Phase boundary extraction ─────────────────────────────────────────────────
Uc_vals = Float64[]
for λ in λ_vals
    uc = 0.0
    for row in kmh_data[λ]
        row.mxy > mxy_threshold && (uc = row.U)
    end
    push!(Uc_vals, uc)
end

λ_all  = λ_vals
Uc_all = Uc_vals

println("Phase boundary Uc vs λ:")
for (λ, uc) in zip(λ_all, Uc_all)
    @printf("  λ = %.2f   Uc = %.3f\n", λ, uc)
end

fig6 = Figure(size=(500, 420))
ax6  = Axis(fig6[1, 1];
    xlabel       = "λ / t",
    ylabel       = "Ũc / t",
    title        = "KMH phase boundary (Rachel & Le Hur Fig. 6)",
    xgridvisible = false,
    ygridvisible = false,
    limits       = ((-0.05, 1.05), (-0.3, 11.0)))

text!(ax6, 0.5, 6.5; text="SDW",   fontsize=20, align=(:center, :center))
text!(ax6, 0.5, 1.5; text="TI/PM", fontsize=20, align=(:center, :center))

scatterlines!(ax6, λ_all, Uc_all;
    color      = :crimson,
    markersize = 10,
    linewidth  = 2,
    overdraw   = true)

display(fig6)
save(joinpath(@__DIR__, "../../docs/src/fig/kmh_phase_diagram.png"), fig6)
println("Saved: docs/src/fig/kmh_phase_diagram.png")
