using VIDA
import Comrade as CM
import CairoMakie as cMakie
using Pyehtim
using VLBISkyModels
import OptimizationBBO as OBBO

file_name = "sgra_Ma+0.5_30_1_0230_1_40_n001.h5"#I0_a_0.94_i_45_sk_1.0_bphi_1.0_br_1.0_bvapp_0_gfact_3_mp_1.833_gp_-1.815_sp_3.574.fits"
in_model = joinpath("data", "selected_sgra_images", file_name)
file = joinpath(dirname(@__DIR__), in_model)
in_img = ehtim.image.load_image(file)
in_img.display()
in_img.blur_circ(μas2rad(20)).display()

sze = pyconvert(Float64, in_img.imvec.size) |> sqrt |> Int
im = reshape(pyconvert(Vector{Float64}, in_img.imvec), (sze, sze))
temp = reshape(map(x->x > sum(im)/(sze^2) ? 1.0 : 0.0, im), (sze, sze))
totx = 0
toty = 0
centervals = []

for x in range(1, sze)
	for y in range(1, sze)
		totx += temp[x, y]*x
		toty += temp[x, y]*y
	end
end
cx = Int(floor(totx / sum(temp)))
cy = Int(floor(toty / sum(temp)))

rat = Int(sze-cy)//Int(sze-cx)
for x in range(1, sze)
	Δx = x-cx
	Δy = sze-cy
    count = 0
	avex = 0
	avey = 0

	for locy in range(cy, sze)
		locx = Int(floor((Δx/Δy)*(locy - cy) + cx))
		curr = im[locx, locy]
		avex += curr*locx
		avey += curr*locy
        count += curr
	end
	append!(centervals, Int.(floor.((avex, avey) ./ count)))

	Δy = 1-cy
    count = 0
    avex = 0
    avey = 0

	for locy in range(1, cy)
		locx = Int(floor((Δx/Δy)*(locy - cy) + cx))
		curr = im[locx, locy]
		avex += curr*locx
		avey += curr*locy
        count += curr
	end
	append!(centervals, Int.(floor.((avex, avey) ./ count)))
end
for y in range(1, sze)
	Δx = sze-cx
	Δy = y-cy
    avex = 0
    avey = 0
    count = 0

	for locx in range(cx, sze)
		locy = Int(floor((Δy/Δx)*(locx - cx) + cy))
		curr = im[locx, locy]
		avex += curr*locx
		avey += curr*locy
        count += curr

	end
	append!(centervals, Int.(floor.((avex, avey) ./ count)))

	Δx = 1-cx
    count = 0
    avex = 0
    avey = 0

	for locx in range(1, cx)
		locy = Int(floor((Δy/Δx)*(locx - cx) + cy))
		curr = im[locx, locy]
		avex += curr*locx
		avey += curr*locy
        count += curr
	end
	append!(centervals, Int.(floor.((avex, avey) ./ count)))

end


begin
	fig = cMakie.Figure(resolution = (800, 800))
	ax = cMakie.Axis(fig[1, 1], titlesize = 20, aspect = 1)
	hm = cMakie.heatmap!(ax, im, colormap = :afmhot)

	hm2 = cMakie.heatmap!(ax, temp, colormap = :afmhot, alpha = 0.3)

	cMakie.scatter!(ax, centervals[1:10:end], centervals[2:10:end], color = :red, markersize = 5)
	cMakie.scatter!(ax, cx, cy, color = :white, markersize = 20, marker = :star5)
	display(fig)
end

