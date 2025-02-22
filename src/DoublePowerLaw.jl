struct DoublePowerLaw{T} <: CM.ComradeBase.AbstractModel
    p1::T
    p2::T
end
CM.ComradeBase.isprimitive(::DoublePowerLaw) = CM.ComradeBase.Isprimitive()
function CM.ComradeBase.intensity_point(m::DoublePowerLaw, p)
    (; p1, p2) = m
    X, Y = p
    r = hypot(X, Y)
    return (r^(-p1)) / (1 + r^(-(p1 + p2)))
end

Base.@constprop :aggressive @inline ComradeBase.visanalytic(
    ::Type{<:DoublePowerLaw{T}},
) where {T} = CM.ComradeBase.NotAnalytic()
Base.@constprop :aggressive @inline ComradeBase.imanalytic(
    ::Type{<:DoublePowerLaw{T}},
) where {T} = CM.ComradeBase.IsAnalytic()
Base.@constprop :aggressive @inline ComradeBase.ispolarized(
    ::Type{<:DoublePowerLaw{T}},
) where {T} = CM.ComradeBase.NotPolarized()
function CM.ComradeBase.radialextent(m::DoublePowerLaw)
    return 2
end
