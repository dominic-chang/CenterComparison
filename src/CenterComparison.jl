using VIDA
import Comrade as CM
using Plots

in_model = "I1_regrid"
file = joinpath((@__DIR__), "..", "data", in_model*".fits")
img = VIDA.load_image(file)

Plots.plot(img)

bh = VIDA.Bhattacharyya(img)


m = CM.MRing(1.0,1.0) # Create ring of 1 radian radius
Plots.plot(m) # plot the ring
m = CM.stretched(m, CM.μas2rad(20), CM.μas2rad(20)) 
m |> Plots.plot #stretch the ring (adjust shape)
m = CM.smoothed(m, CM.μas2rad(10)/(2*√(2*log(2)))) 
m |> Plots.plot #smooth the ring (adjust blur)
m = CM.shifted(m, μas2rad(5), μas2rad(2)) 
m |> Plots.plot #shift the ring (adjust the position)


# Create Our Ring Model
# Params: rad, width, α, β, ell, x, y

function ring_model(params)
    (;rad, width, α, β, ell, x, y) = params
    m = CM.MRing(α, β) # Create ring of 1 radian radius
    m = CM.stretched(m, CM.μas2rad(rad), CM.μas2rad(rad*ell)) 
    m = CM.smoothed(m, CM.μas2rad(width)/(2*√(2*log(2)))) 
    return CM.ThreadedModel(CM.shifted(m, μas2rad(x), μas2rad(y)))
end

VIDA.divergence(bh, m)

#lower bounds on parameters
lower = (
    rad = (8),
    width = (5),
    α = -1,
    β = -1,
    ell = 0,
    x = (-10),
    y = (-10),
)

#upper bounds on parameters
upper = (
    rad = (40),
    width = (40),
    α = 1,
    β = 1,
    ell = 1,
    x = (10),
    y = (10),
)
# Defines a VIDAproblem type for optimization
# Takes a divergenced function (bh) and a model function (ring_model)
# Our choice of optimization requires a search range with lower and upper bounds for the ring model parameters
prob = VIDA.VIDAProblem(bh, ring_model, lower, upper) 
f, t, (lb, ub) = VIDA.build_opt(prob, true)


import OptimizationBBO as OBBO
count = 0
function callback(state, loss_val; doplot = false) # callback function to print iterations and loss
    global count
    count += 1
    if count %100 |> iszero
        println("$count: $loss_val")
    end
    return false
end
xopt, opt_temp, divmin = vida(prob, OBBO.BBO_adaptive_de_rand_1_bin(); maxiters=50_000, callback=callback) #run problem 
gr = VIDA.imagepixels(μas2rad(100), μas2rad(100), 100,100)
VIDA.intensitymap(opt_temp, gr) |> Plots.plot
xopt
fig = triptic(img, opt_temp)

# Save results
outpath = joinpath((@__DIR__),"..","results")
try
    mkdir(outpath)
catch e
    println(e)
end
fileout = open(joinpath(outpath, in_model), "w")

write(fileout, "upper = " * string(upper) * "\n")
write(fileout, "lower = " * string(lower) * "\n")
write(fileout, "best_fit = " * string(xopt))
close(fileout)
Plots.savefig(fig, joinpath(outpath, in_model*"_triptic.png"))



