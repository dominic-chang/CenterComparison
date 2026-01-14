using VIDA
import Comrade as CM
import CairoMakie as cMakie
using Pyehtim
using VLBISkyModels
import OptimizationBBO as OBBO

file_name = "sgra_Ma+0.5_30_1_0230_1_40_n000.h5"#I0_a_0.94_i_45_sk_1.0_bphi_1.0_br_1.0_bvapp_0_gfact_3_mp_1.833_gp_-1.815_sp_3.574.fits"
in_model = joinpath("data", "selected_sgra_images", file_name)
file = joinpath(dirname(@__DIR__), in_model)
in_img = ehtim.image.load_image(file)
in_img.display()
in_img.blur_circ(μas2rad(20)).display()

using Krang
metric = Krang.Kerr(0.94);

θo = 17 * π / 180;
ρmax = 12.0;
camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, 400);

χ = -1.7;
ι = 0.58;
βv = 0.87;
σ = 0.73;
η1 = 2.64;
η2 = π - η1;
n = 1

magfield1 = Krang.SVector(sin(ι) * cos(η1), sin(ι) * sin(η1), cos(ι));
vel = Krang.SVector(βv, (π / 2), χ);
R = 4.0;
p1 = 4.0;
p2 = 4.0;

material1 = Krang.ElectronSynchrotronPowerLawPolarization(
	magfield1...,
	vel...,
	σ,
	R,
	p1,
	p2,
	(n,),
);
θs = (90 * π / 180);
geometry1 = Krang.ConeGeometry(θs)

mesh1 = Krang.Mesh(geometry1, material1)

scene = Krang.Scene((mesh1, ))
stokesvals = render(camera, scene)

im = map(x->x[1], stokesvals)
sze = size(im)[1]
#sze = pyconvert(Float64, in_img.imvec.size) |> sqrt |> Int
#im = reshape(pyconvert(Vector{Float64}, in_img.imvec), (sze, sze))
temp = reshape(map(x->x > sum(im)/(sze^2) ? 1.0 : 0.0, im), (sze, sze))
begin
	totx = 0
	toty = 0

	for x in range(1, sze)
		for y in range(1, sze)
			totx += temp[x, y]*x
			toty += temp[x, y]*y
		end
	end
end
cx = Int(floor(totx / sum(temp)))
cy = Int(floor(toty / sum(temp)))
centervals = []


rat = Int(sze-cy)//Int(sze-cx)
for x in range(1, sze)
	Δx = x-cx
	Δy = sze-cy
	avex = 0
	avey = 0
    count = 0

	for locy in range(cy, sze)
		locx = max(Int(floor((Δx/Δy)*(locy - cy) + cx)), 1)
		if temp[locx, locy] == 1
			avex += locx
			avey += locy
            count += 1
		end
	end
    count > 0 && append!(centervals, Int.(floor.((avex, avey) ./ count)))

	Δy = 1-cy
    count = 0
    avex = 0
    avey = 0
	for locy in range(1, cy)
		locx = max(Int(floor((Δx/Δy)*(locy - cy) + cx)), 1)
		if temp[locx, locy] == 1
			avex += locx
			avey += locy
            count += 1
		end
	end
    count > 0 && append!(centervals, Int.(floor.((avex, avey) ./ count)))
end
for y in range(1, sze)
	Δx = sze-cx
	Δy = y-cy
    avex = 0
    avey = 0
    count = 0

	for locx in range(cx, sze)
		locy = max(Int(floor((Δy/Δx)*(locx - cx) + cy)), 1)
		if temp[locx, locy] == 1
			avex += locx
			avey += locy
            count += 1
		end
	end
    count > 0 && append!(centervals, Int.(floor.((avex, avey) ./ count)))


	Δx = 1-cx
    count = 0
    avex = 0
    avey = 0

	for locx in range(1, cx)
		locy = max(Int(floor((Δy/Δx)*(locx - cx) + cy)), 1)
		if temp[locx, locy] == 1
			avex += locx
			avey += locy
            count += 1
		end
	end
    count > 0 && append!(centervals, Int.(floor.((avex, avey) ./ count)))


end


radii = map(i->emission_radius(camera.screen.pixels[centervals[2i-1], centervals[2i]], θs, true, n)[1] + emission_radius(camera.screen.pixels[centervals[2i-1], centervals[2i]], θs, false, n)[1], range(1, length(centervals)÷2))
begin
	fig = cMakie.Figure(resolution = (800, 800))
	ax = cMakie.Axis(fig[1, 1], titlesize = 20, aspect = 1)
	hm = cMakie.heatmap!(ax, im, colormap = :afmhot)

	hm2 = cMakie.heatmap!(ax, temp, colormap = :afmhot, alpha = 0.3)

	cMakie.scatter!(ax, centervals[1:10:end], centervals[2:10:end], color = :red, markersize = 5)
	cMakie.scatter!(ax, cx, cy, color = :white, markersize = 20, marker = :star5)
    cMakie.scatter!(ax, GI.x(cent), GI.y(cent), color = :red) 
	display(fig)
end

cMakie.scatter(radii)