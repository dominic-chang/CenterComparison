using VIDA
import Comrade as CM
import CairoMakie as cMakie
using Pyehtim
using VLBISkyModels
import OptimizationBBO as OBBO

for file_name in readdir(joinpath(dirname(@__DIR__), "thickRings"))
    in_model = joinpath("thickRings", file_name)
    file = joinpath(dirname(@__DIR__), in_model)
    in_img = ehtim.image.load_image(file)
    in_img.display()

    sze = pyconvert(Float64, in_img.imvec.size) |> sqrt |> Int
    fovx = pyconvert(Float64, in_img.fovx())
    fovy = pyconvert(Float64, in_img.fovy())
    img = IntensityMap(
        reverse(reshape(pyconvert(Vector{Float64}, in_img.imvec), (sze, sze)), dims = 1),
        CM.imagepixels(fovx, fovy, sze, sze),
    )

    cMakie.plot(img)

    grid = VIDA.imagepixels(μas2rad(100), μas2rad(100), 100, 100)
    m = CM.Crescent(1.0, 0.5, 0.0, 0.0) # Create ring of 1 radian radius
    m = CM.stretched(m, CM.μas2rad(20), CM.μas2rad(20))
    cMakie.image(intensitymap(m, grid)) # plot the ring
    m = CM.shifted(m, μas2rad(5), μas2rad(2))
    cMakie.image(intensitymap(m, grid)) # plot the ring


    # Create Our Ring Model
    # Params: rad, width, α1, β1, α2, β2, ell, x, y

    function ring_model(params)
        (; rad, rad_rat, shift, rot, ell, x, y) = params
        m = CM.Crescent(1.0, rad_rat, shift * (1 - rad_rat), 0.0) # Create ring of 1 radian radius
        #m = CM.stretched(m, CM.μas2rad(rad), CM.μas2rad(rad * ell))
        #return CM.shifted(m, μas2rad(x), μas2rad(y))
        return CM.modify(
            m,
            Rotate(rot),
            Stretch(CM.μas2rad(rad), CM.μas2rad(rad * ell)),
            Shift(μas2rad(x), μas2rad(y)),
        )
    end
    #lower bounds on parameters
    lower_bound =
        (rad = 8.0, rad_rat = 0.1, shift = 0.0, rot = 0.0, ell = 0.0, x = -10.0, y = -10.0)

    #upper bounds on parameters
    upper_bound =
        (rad = 50.0, rad_rat = 0.99, shift = 1.0, rot = 2π, ell = 1.0, x = 10.0, y = 10.0)
    bh = VIDA.NxCorr(img) # Loss function for comparison

    # Defines a VIDAproblem type for optimization
    # Takes a divergenced function (bh) and a model function (ring_model)
    # Our choice of optimization requires a search range with lower and upper bounds for the ring model parameters
    prob = VIDA.VIDAProblem(bh, ring_model, lower_bound, upper_bound)
    f, t, (lb, ub) = VIDA.build_opt(prob, true)


    global count = 0
    function callback(state, loss_val; doplot = false) # callback function to print iterations and loss
        global count
        count += 1
        if count % 100 |> iszero
            println("$count: $loss_val")
        end
        return false
    end
    xopt, opt_temp, divmin = vida(
        prob,
        OBBO.BBO_adaptive_de_rand_1_bin();
        maxiters = 50_000,
        callback = callback,
    ) #run problem 
    gr = VIDA.imagepixels(μas2rad(100), μas2rad(100), 100, 100)
    img_temp = VIDA.intensitymap(opt_temp, gr)
    img_temp |> cMakie.plot
    xopt
    fig = triptic(img, opt_temp)

    # Save results
    outpath = joinpath((@__DIR__), "..", "results", "crescent")
    mkpath(outpath)
    fileout = open(joinpath(outpath, file_name .* ".txt"), "w")

    write(fileout, "upper = " * string(upper_bound) * "\n")
    write(fileout, "lower = " * string(lower_bound) * "\n")
    write(fileout, "best_fit = " * string(xopt))
    close(fileout)
    cMakie.save(joinpath(outpath, file_name * "_triptic.png"), fig)

end
