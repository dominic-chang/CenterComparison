using VIDA

include(joinpath((@__DIR__), "generate_n_ring.jl"))

θo = 30 / 180 * π
a = 0.94
rs = 4.0

n = 1
fig = Figure();
ax = Axis(fig[1, 1], aspect = 1.0)
xmin = -10
xmax = 10
ymin = -10
ymax = 10
xlims!(ax, xmin, xmax)
ylims!(ax, ymin, ymax)

bulkres = 400
screenres = 200

mat = generate_n_ring(
    rs,
    θo,
    a,
    n,
    bulkres,
    screenres,
    xmin,
    xmax,
    ymin,
    ymax;
    tolerance = 1e-6,
    alpha = 1e-5,
)
αβvals = get_isoradial_curve(rs, θo, a, n, bulkres; tolerance = 1e-6)
ρvals = [hypot(x...) for x in αβvals]
φvals = [atan(x[2], x[1]) for x in αβvals]
fig = Figure()
ax = Axis(fig[1, 1], aspect = 1)
rsvals = []
for i = 1:length(ρvals)
    p = [φvals[i], θo, a, n]
    x = ρvals[i]
    append!(rsvals, rad(x, p))
end
scatter(rsvals ./ rs)

ax = Axis(fig[1, 1], aspect = 1)
heatmap!(
    ax,
    range(xmin, xmax, length = screenres),
    range(ymin, ymax, length = screenres),
    mat',
)
scatter!(ax, αβvals)


display(fig)

grid = imagepixels(2xmax * μas2rad(5.03), 2ymax * μas2rad(5.03), screenres, screenres)
intmap = IntensityMap(mat', grid)
save_fits(
    joinpath(
        dirname((@__DIR__)),
        "thinRings",
        "m_d_5_03_rs_$(rs)_inc_$(Int(round(θo*180/pi)))_a_$(round(a;digits=2))_$(n)_ring.fits",
    ),
    intmap,
)
intmap |> imageviz
