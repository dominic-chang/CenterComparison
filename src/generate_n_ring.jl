using Krang
using Optimization, OptimizationOptimJL
using Enzyme
using CairoMakie

Δ(r, a) = r^2 - 2r + a^2

function rad(ρ, p)
    φ, θo, a, n = p
    β, α = ρ .* sincos(φ)
    met = Krang.Kerr(a)
    pix = Krang.IntensityPixel(met, α, β, θo)
    Krang.emission_radius(pix, π / 2, false, n)[1] +
    Krang.emission_radius(pix, π / 2, true, n)[1]
end

function rmp(a)
    return 2 .* (1 .+ cos.(2 / 3 .* acos.(a .* (-1, 1))))
end

function λcrit(r, a)
    return a + r / a * (r - 2 * Δ(r, a) / (r - 1))
end

function ηcrit(r, a)
    return r^3 / a^2 * (4 * Δ(r, a) / (r - 1)^2 - r)
end

function ληcrit(r, a)
    return (
        a + r / a * (r - 2 * Δ(r, a) / (r - 1)),
        r^3 / a^2 * (4 * Δ(r, a) / (r - 1)^2 - r),
    )
end

function αβ(λη, a, θo)
    λ, η = λη
    return (-λ / sin(θo), √(η + a^2 * cos(θo)^2 - λ^2 * cot(θo)^2))
end

function get_isoradial_curve(θo, a, res)
    #generate critical values to sample from
    rcritvals = range(rmp(a)..., length = res)
    ληcritvals =
        filter(x -> x[2] > -(a^2 * cos(θo)^2 - x[1]^2 * cot(θo)^2), ληcrit.(rcritvals, a))
    αβcritvals1 = αβ.(ληcritvals, a, θo)
    αβcritvals2 = map(x -> (x[1], -x[2]), αβcritvals1)
    αβcritvals = append!(αβcritvals1, αβcritvals2[end:-1:begin])

    ρφcritvals = map(x -> (hypot(x...), atan(x[end:-1:begin]...)), (αβcritvals))
    return ρφcritvals
end

function get_isoradial_curve(rs, θo, a, n, res; tolerance = 1.001)
    ρφcritvals = get_isoradial_curve(θo, a, res)
    ρφvals = map(
        x -> begin
            ρ, φ = x
            p = [φ, θo, a, n]
            x0 = [tolerance * ρ]
            optprob = OptimizationFunction(
                (x, p) -> (rs - rad(x, p))^2,
                Optimization.AutoEnzyme(),
            )
            prob = Optimization.OptimizationProblem(optprob, x0, p)
            ans = solve(prob, Adam(; alpha = 1e-5))
            (ans, φ)
        end,
        ρφcritvals,
    )
    return map(x -> (x[1][1] * cos(x[2]), x[1][1] * sin(x[2])), ρφvals)
end

function generate_n_ring(rs, θo, a, n, bulkres, screenres, xmin, xmax, ymin, ymax)
    αβvals = get_isoradial_curve(rs, θo, a, n, bulkres)
    toppoints = sort(filter(x -> x[2] > 0, αβvals), by = x -> x[1])
    bottompoints = sort(filter(x -> x[2] <= 0, αβvals), by = x -> x[1])

    toppoints = sort(filter(x -> x[2] > 0, αβvals), by = x -> x[1])
    bottompoints = sort(filter(x -> x[2] <= 0, αβvals), by = x -> x[1])
    leftedgepoints = [toppoints[1], bottompoints[1]]
    rightedgepoints = [bottompoints[end], toppoints[end]]

    mat = zeros(screenres, screenres)

    for points in (toppoints, leftedgepoints, bottompoints, rightedgepoints)
        currpoint, nextpoint = points[1:2]
        m = (nextpoint[2] - currpoint[2]) / (nextpoint[1] - currpoint[1])
        yo = currpoint[2] - m * currpoint[1]
        xo = currpoint[1] - currpoint[2] / m
        Δα = (xmax - xmin) / (screenres - 1)
        Δβ = (ymax - ymin) / (screenres - 1)
        curr = 2

        done = false
        col = 0
        while col < screenres && !done
            col += 1
            row = 0
            while row < screenres && !done
                row += 1

                αmin = (col - 1) * Δα + xmin
                αmax = Δα + αmin

                βmin = (row - 1) * Δβ + ymin
                βmax = Δβ + βmin

                p1, p2 = sort([currpoint, nextpoint], by = x -> x[1])
                βb, βt = minmax(p1[2], p2[2])
                if αmin > p2[1]
                    curr += 1
                    col -= 1
                    row = 0
                    if curr > length(points)
                        done = true
                        break
                    end
                    currpoint = nextpoint
                    nextpoint = points[curr]
                    m = (nextpoint[2] - currpoint[2]) / (nextpoint[1] - currpoint[1])
                    yo = currpoint[2] - m * currpoint[1]
                    xo = currpoint[1] - currpoint[2] / m
                    p1, p2 = sort([currpoint, nextpoint], by = x -> x[1])
                    βb, βt = minmax(p1[2], p2[2])
                end

                # Check to see if line between points intersects cell
                if (
                    (βmin <= (m * αmin + yo) <= βmax) ||
                    (βmin <= (m * αmax + yo) <= βmax) ||
                    (αmin <= (βmin / m + xo) <= αmax) ||
                    (αmin <= (βmax / m + xo) <= αmax)
                ) && !(αmax < p1[1] || αmin > p2[1] || βmin > βt || βmax < βb)
                    mat[row, col] = 1
                end
            end
        end
    end
    return mat
end
