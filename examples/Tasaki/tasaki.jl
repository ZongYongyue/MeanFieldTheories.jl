using MeanFieldTheories
using CairoMakie
using LinearAlgebra, Printf

t = 1.0
ds = [1.0, 2.0/sqrt(3), 1.5]
Ua = t
Ub = t
mua = 2t

a1 = [1.0,  0.0]
a2 = [0.5, √3/2]

unitcell = Lattice(
    [Dof(:cell, 1), Dof(:sub, 2, [:A, :B])],
    [QN(cell=1, sub=1), QN(cell=1, sub=2)],
    [[0.0, 0.0], [0.5, √3/2]];
    vectors = [a1, a2]
)

dofs = SystemDofs([Dof(:cell, 1), Dof(:sub, 2, [:A, :B]), Dof(:spin, 2, [:up, :dn])])

nn_bonds     = bonds(unitcell, (:p, :o), 1)
onsite_bonds = bonds(unitcell, (:p, :o), 0)

twobody = generate_twobody(dofs, onsite_bonds,
    (deltas, qn1, qn2, qn3, qn4) -> begin
        if qn1.sub == qn2.sub == 1
            if (qn1.spin == qn2.spin) && (qn3.spin == qn4.spin) && (qn1.spin !== qn3.spin)
                return Ua/2
            else
                return 0.0
            end
        elseif qn1.sub == qn2.sub == 2
            if (qn1.spin == qn2.spin) && (qn3.spin == qn4.spin) && (qn1.spin !== qn3.spin)
                return Ub/2
            else
                return 0.0
            end
        else
            return 0.0
        end
    end,
    order = (cdag, :i, c, :i, cdag, :i, c, :i))

kpoints = [[r, 0.0] for r in range(0, 2pi, length=65)[1:end-1]]
n_elec = length(kpoints)

labels = ["d = 1.0", "d = 2/√3", "d = 1.5"]

fig = Figure(size=(1500, 500))
kx = collect(range(0, 2pi, length=65)[1:end-1])

for (id, d) in enumerate(ds)
    mub = t * d^2

    onebody_h = generate_onebody(dofs, nn_bonds,
        (delta, qn1, qn2) -> begin
            if qn1.spin !== qn2.spin
                return 0.0
            else
                if qn1.sub == qn2.sub == 2
                    return 0.0
                elseif qn1.sub == qn2.sub == 1
                    return t
                elseif qn1.sub !== qn2.sub
                    return t*d
                end
            end
        end)

    onebody_m = generate_onebody(dofs, onsite_bonds,
        (delta, qn1, qn2) -> begin
            if qn1.spin !== qn2.spin
                return 0.0
            else
                if qn1.sub == qn2.sub == 1
                    mua
                elseif qn1.sub == qn2.sub == 2
                    mub
                else
                    throw(ArgumentError("error"))
                end
            end
        end,
    order = (cdag, :i, c, :i))

    onebody = (ops=[onebody_h.ops; onebody_m.ops], irvec=[onebody_h.irvec; onebody_m.irvec], delta=[onebody_h.delta; onebody_m.delta])

    println("\n" * "="^60)
    println("  d = $d  ($(labels[id]))")
    println("="^60)

    hf = solve_hfk(dofs, onebody, twobody, kpoints, n_elec;
        n_restarts     = 10,
        G_init         = nothing,
        field_strength = 1.0,
        n_warmup       = 15,
        tol            = 1e-12,
        verbose        = true)

    mags = local_spin(dofs, hf.G_k)
    print_spin(mags)

    ph = solve_ph_excitations(dofs, onebody, twobody, hf, kpoints, [[pi*2.0, 0.0]];
        n_list  = [1,2],
        solver = :TDA,
        verbose = true)

    ax = Axis(fig[1, id]; xlabel="q", ylabel="Energy", title=labels[id])
    n_bands = length(ph.energies[1])
    for band in 1:n_bands
        E_band = [sort(ph.energies[qi])[band] for qi in 1:length(kpoints)]
        scatter!(ax, kx, E_band; color=:blue, markersize=3)
    end
    ylims!(ax, 0, 0.6)
end

display(fig)

save("docs/src/fig/tasaki.png", fig)
