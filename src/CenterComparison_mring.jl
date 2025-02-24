using VIDA
import Comrade as CM
using Plots
import CairoMakie as cMakie
using Pyehtim

for file_name in
    readdir(joinpath(dirname(@__DIR__), "data", "selected_sgra_images"))[185:end]
    in_model = joinpath("selected_sgra_images", file_name)
    file = joinpath((@__DIR__), "..", "data", in_model)
    in_img = ehtim.image.load_image(file)
    in_img.display()

    sze = pyconvert(Float64, in_img.imvec.size) |> sqrt |> Int
    fovx = pyconvert(Float64, in_img.fovx())
    fovy = pyconvert(Float64, in_img.fovy())
    img = IntensityMap(
        reverse(reshape(pyconvert.(Float64, in_img.imvec), (sze, sze)), dims = 1),
        CM.imagepixels(fovx, fovy, sze, sze),
    )

    Plots.plot(img)

    m = CM.MRing(1.0, 1.0) # Create ring of 1 radian radius
    Plots.plot(m) # plot the ring
    m = CM.stretched(m, CM.μas2rad(20), CM.μas2rad(20))
    m |> Plots.plot #stretch the ring (adjust shape)
    m = CM.smoothed(m, CM.μas2rad(10) / (2 * √(2 * log(2))))
    m |> Plots.plot #smooth the ring (adjust blur)
    m = CM.shifted(m, μas2rad(5), μas2rad(2))
    m |> Plots.plot #shift the ring (adjust the position)


    # Create Our Ring Model
    # Params: rad, width, α1, β1, α2, β2, ell, x, y

    function ring_model(params)
        (; rad, width, α1, β1, α2, β2, ell, x, y) = params
        m = CM.MRing([α1, α2], [β1, β2]) # Create ring of 1 radian radius
        m = CM.stretched(m, CM.μas2rad(rad), CM.μas2rad(rad * ell))
        m = CM.smoothed(m, CM.μas2rad(width) / (2 * √(2 * log(2))))
        return CM.shifted(m, μas2rad(x), μas2rad(y))
    end
    #lower bounds on parameters
    lower_bound = (
        rad = 8.0,
        width = 0.5,
        α1 = -1.0,
        β1 = -1.0,
        α2 = -1.0,
        β2 = -1.0,
        ell = 0.0,
        x = -10.0,
        y = -10.0,
    )

    #upper bounds on parameters
    upper_bound = (
        rad = 50.0,
        width = 80.0,
        α1 = 1.0,
        β1 = 1.0,
        α2 = 1.0,
        β2 = 1.0,
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
    img_temp =VIDA.intensitymap(opt_temp, gr) 
    img_temp |> Plots.plot
    xopt
    fig = triptic(img, opt_temp)

    # Save results
    outpath = joinpath((@__DIR__), "..", "results")
    try
        mkdir(outpath)
    catch e
        println(e)
    end
    fileout = open(joinpath(outpath, "mring", in_model .* ".txt"), "w")

    write(fileout, "upper = " * string(upper_bound) * "\n")
    write(fileout, "lower = " * string(lower_bound) * "\n")
    write(fileout, "best_fit = " * string(xopt))
    close(fileout)
    cMakie.save(joinpath(outpath, in_model * "_triptic.png"), fig)

end
