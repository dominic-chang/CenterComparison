using Krang

struct UniformEmission{T} <: Krang.AbstractMaterial
    subimgs::Tuple{Int}
    rmin::T
    rmax::T
end

function (m::UniformEmission)(
    pix::Krang.AbstractPixel{T},
    intersection::Krang.Intersection,
) where {T}
    (; rs) = intersection

    return m.rmin < rs < m.rmax
end

a = 0.001
metric = Krang.Kerr(a);
θo = 30 * π / 180;
ρmax = 10.0;
rmin = 4.0; # minimum radius to truncate cone
rmax = 5.0; # maximum radius to truncate cone
n = 0; # sub-image to ray trace
screenres = 200

camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, screenres);
mesh = Krang.Mesh(Krang.ConeGeometry((π / 2.0)), UniformEmission((n,), rmin, rmax))
scene = Krang.Scene((mesh,))

mat = render(camera, scene)

import CairoMakie as CMk

theme = CMk.Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
    ),
)

CMk.set_theme!(CMk.merge!(theme, CMk.theme_latexfonts()))

fig = CMk.Figure(resolution = (700, 700));
ax = CMk.Axis(fig[1, 1], titlesize = 20, aspect = 1)
hm = CMk.heatmap!(ax, mat, colormap = :afmhot)
CMk.Colorbar(fig[1, 2], hm, label = "Intensity", labelsize = 20)
display(fig)

using VIDA
grid = imagepixels(2ρmax * μas2rad(5.03), 2ρmax * μas2rad(5.03), screenres, screenres)
intmap = IntensityMap(mat, grid)
save_fits(
    joinpath(
        dirname((@__DIR__)),
        "thickRings",
        "m_d_5_03_rs_$(rmin)_$(rmax)_inc_$(Int(round(θo*180/pi)))_a_$(round(a;digits=2))_$(n)_ring.fits",
    ),
    intmap,
)
intmap |> imageviz
