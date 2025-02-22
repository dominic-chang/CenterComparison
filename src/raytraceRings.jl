include(joinpath((@__DIR__), "generate_n_ring.jl"))

θo = 89 / 180 * π
a = 0.999
rs = 4.2

n = 10
fig = Figure();
ax = Axis(fig[1, 1], aspect = 1.0)
xmin = -10
xmax = 10
ymin = -10
ymax = 10
xlims!(ax, xmin, xmax)
ylims!(ax, ymin, ymax)

bulkres = 500
screenres = 300

mat = generate_n_ring(rs, θo, a, n, bulkres, screenres, xmin, xmax, ymin, ymax)
αβvals = get_isoradial_curve(rs, θo, a, n, bulkres)

ax = Axis(fig[1, 1], aspect = 1)
heatmap!(
    ax,
    range(xmin, xmax, length = screenres),
    range(ymin, ymax, length = screenres),
    mat',
)
scatter!(ax, αβvals)
display(fig)
