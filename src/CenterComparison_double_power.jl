using VIDA
import Comrade as CM
using Plots
include(joinpath((@__DIR__), "DoublePowerLaw.jl"))

in_model = "I1_a_0.94_i_45_sk_1.0_bphi_1.0_br_1.0_bvapp_0_gfact_3_mp_1.833_gp_-1.815_sp_3.574"
file = joinpath((@__DIR__), "..", "data", in_model*".fits")
img = VIDA.load_image(file)

Plots.plot(img)

m = DoublePowerLaw(1.0,5.0) # Create ring of 1 radian radius
Plots.plot(m) # plot the ring
m = CM.stretched(m, CM.μas2rad(20), CM.μas2rad(20)) 
m |> Plots.plot #stretch the ring (adjust shape)
m = CM.shifted(m, μas2rad(5), μas2rad(2)) 
m |> Plots.plot #shift the ring (adjust the position)

function ring_model(params)
    (;rad, p1, p2, ell, x, y) = params
    m = DoublePowerLaw(p1, p2) # Create ring of 1 radian radius
    m = CM.stretched(m, CM.μas2rad(rad), CM.μas2rad(rad*ell)) 
    return CM.ThreadedModel(CM.shifted(m, μas2rad(x), μas2rad(y)))
end

#lower bounds on parameters
lower_bound = (
    rad = 8.0,
    p1 = 1.0,
    p2 = 0.1,
    ell = 0.0,
    x = -10.0,
    y = -10.0,
)

#upper bounds on parameters
upper_bound = (
    rad = 50.0,
    p1 = 100.0,
    p2 = 100.0,
    ell = 1.0,
    x = 10.0,
    y = 10.0,
)
bh = VIDA.NxCorr(img) # Loss function for comparison

# Defines a VIDAproblem type for optimization
# Takes a divergenced function (bh) and a model function (ring_model)
# Our choice of optimization requires a search range with lower and upper bounds for the ring model parameters
prob = VIDA.VIDAProblem(bh, ring_model, lower_bound, upper_bound) 
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
xopt, opt_temp, divmin = vida(prob, OBBO.BBO_adaptive_de_rand_1_bin(); maxiters=100_000, callback=callback) #run problem 
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
fileout = open(joinpath(outpath, in_model.*".txt"), "w")

write(fileout, "upper = " * string(upper_bound) * "\n")
write(fileout, "lower = " * string(lower_bound) * "\n")
write(fileout, "best_fit = " * string(xopt))
close(fileout)
Plots.savefig(fig, joinpath(outpath, in_model*"_triptic.png"))



