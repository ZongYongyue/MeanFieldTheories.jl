"""
Plot HF band structure and magnon dispersion for the t-U-J-J' model.

Run:
    julia --project=examples examples/Heisenberg/plot.jl
"""

using CairoMakie
using Printf

# ── Helper: parse a .dat file ─────────────────────────────────────────────────
function read_dat(path)
    lines_raw = readlines(path)
    xtick_pos = Float64[]
    xtick_lab = String[]
    mu        = nothing
    for line in lines_raw
        if startswith(line, "# xtick_pos:")
            s = strip(split(line, "xtick_pos:")[2])
            xtick_pos = parse.(Float64, split(strip(s, ['[', ']']), ','))
        elseif startswith(line, "# xtick_lab:")
            s  = strip(split(line, "xtick_lab:")[2])
            s2 = replace(s, "[" => "", "]" => "", "\"" => "")
            xtick_lab = strip.(split(s2, ','))
        elseif startswith(line, "# mu:")
            mu = parse(Float64, strip(split(line, "mu:")[2]))
        end
    end
    data_lines = filter(l -> !startswith(l, "#") && !isempty(strip(l)), lines_raw)
    arc_vals  = Float64[]
    all_bands = Vector{Float64}[]
    for line in data_lines
        cols = split(strip(line))
        push!(arc_vals, parse(Float64, cols[1]))
        push!(all_bands, parse.(Float64, cols[2:end]))
    end
    return arc_vals, all_bands, xtick_pos, xtick_lab, mu
end

# ── Read data files ───────────────────────────────────────────────────────────
hf_file = joinpath(@__DIR__, "hf_bands.dat")
isfile(hf_file) || error("hf_bands.dat not found — run run.jl first")
hf_arc, hf_bands_v, xtick_pos, xtick_lab, mu = read_dat(hf_file)
hf_mat = reduce(hcat, hf_bands_v)'   # (Nq, norb)

magnon_file = joinpath(@__DIR__, "magnon.dat")
isfile(magnon_file) || error("magnon.dat not found — run run.jl first")
mag_arc, magnon_bands_v, _, _, _ = read_dat(magnon_file)
mag_mat = reduce(hcat, magnon_bands_v)'   # (Nq, n_modes)

# ── Figure: two panels side by side ───────────────────────────────────────────
fig = Figure(size = (1100, 480))

# ── Panel 1: HF bands ────────────────────────────────────────────────────────
ax1 = Axis(fig[1, 1];
    xlabel       = "",
    ylabel       = "Energy (t)",
    title        = "HF bands (t-U-J-J')",
    xticks       = (xtick_pos, xtick_lab),
    xgridvisible = false,
)
for xp in xtick_pos
    vlines!(ax1, xp; color = :gray60, linewidth = 0.8)
end
for ni in axes(hf_mat, 2)
    lines!(ax1, hf_arc, hf_mat[:, ni]; color = :steelblue, linewidth = 1.5,
           label = ni == 1 ? "HF bands" : nothing)
end
if mu !== nothing
    hlines!(ax1, mu; color = :red, linewidth = 1.2, linestyle = :dash,
            label = @sprintf("μ = %.4f t", mu))
end
axislegend(ax1; position = :rt, merge = true)
xlims!(ax1, hf_arc[1], hf_arc[end])

# ── Panel 2: Magnon dispersion ────────────────────────────────────────────────
ax2 = Axis(fig[1, 2];
    xlabel       = "",
    ylabel       = "Energy (t)",
    title        = "Magnon dispersion (RPA)",
    xticks       = (xtick_pos, xtick_lab),
    xgridvisible = false,
)
for xp in xtick_pos
    vlines!(ax2, xp; color = :gray60, linewidth = 0.8)
end
n_plot = min(2, size(mag_mat, 2))
for ni in 1:n_plot
    lines!(ax2, mag_arc, mag_mat[:, ni]; color = :steelblue, linewidth = 1.2,
           label = ni == 1 ? "magnon" : nothing)
end
axislegend(ax2; position = :rt, merge = true)
xlims!(ax2, mag_arc[1], mag_arc[end])

display(fig)
save(joinpath(@__DIR__, "../../docs/src/fig/heisenberg_magnon.png"), fig; px_per_unit = 2)
println("Saved: docs/src/fig/heisenberg_magnon.png")
