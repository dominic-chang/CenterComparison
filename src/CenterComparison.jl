using VIDA
import Comrade as CM
using Plots
using Pyehtim


file = joinpath((@__DIR__), "..", "data","I1.fits")
img = VIDA.load_image(file)
img |> Plots.plot

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
    rad = μas2rad(8),
    width = μas2rad(5),
    α = -1,
    β = -1,
    ell = 0,
    x = μas2rad(-10),
    y = μas2rad(-10),
)

#upper bounds on parameters
upper = (
    rad = μas2rad(40),
    width = μas2rad(20),
    α = 1,
    β = 1,
    ell = 1,
    x = μas2rad(10),
    y = μas2rad(10),
)
# Defines a VIDAproblem type for optimization
# Takes a divergenced function (bh) and a model function (ring_model)
# Our choice of optimization requires a search range with lower and upper bounds for the ring model parameters
prob = VIDA.VIDAProblem(bh, ring_model, lower, upper) 
f, t, (lb, ub) = VIDA.build_opt(prob, true)


import OptimizationBBO as OBBO
count = 0
function callback(state, loss_val; doplot = false)
    global count
    count += 1
    if count %100 |> iszero
        println("$count: $loss_val")
    end
    return false
end
xopt, opt_temp, divmin = vida(prob, OBBO.BBO_adaptive_de_rand_1_bin(); maxiters=50_000, callback=callback) #run problem 




