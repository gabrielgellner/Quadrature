#==
# Some basic numerical integration routines
> Code taken from: https://github.com/cmcbride/julia/blob/d54095465da8ad06d1d473f6709df0d599f02c7d/base/integrate.jl
Issues to try and fix:
- This assumes that eltype exists, which may not be true for iterables. It also assumes that zero exists, which is not true if e.g. eltype is an array.
- The right thing to do is to consider three cases separately:
    - zero or one elements: return sum(y)*(1.0*dx). The sum function will do the right thing,
      but you need the 1.0*dx to ensure that the type is the same as below. 2 elements: (y[1]+y[2])*(dx*0.5)
    - 3+ elements: iterate, but initialize the summand to y[1]*0.5 + y[2]
      However, you have to use start/next/done rather than length(y) or y[i] to get the length and elements, respectively.

# Future ideas
Port over as much as I can from scipy.integrate. It would be nice to have feature parity
with that package.
==#

# simple integration routines mostly for tabulated data

function trapz(y)
    trapz(y, one(eltype(y)))
end

# integration for uniform points (in x)
function trapz(y, dx::Number)
    r = zero(zero(eltype(y)) * zero(dx))
    r += (y[1] + y[end])/2
    r += sum(y[2:end-1])
    r * dx
end

# integration for non-uniform points (in x)
function trapz(y, x)
    n = length(y)
    if n != length(x)
        throw(ArgumentError("Input x,y must be of same length"))
    end
    r = zero(zero(eltype(x))*zero(eltype(y)))
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0 #TODO: this forces the output to be Float64
end

function simps(y)
    simps(y, one(eltype(y)))
end

function simps(y, dx::Number)
    n = length(y)
    if iseven(n)
        throw(ArgumentError("Simpson rule requires ODD length input (EVEN number of intervals)"))
    end
    r = zero(zero(eltype(y))*zero(dx))
    for i in 2:2:n-1
        r += y[i]
    end
    r *= 2.0
    for i in 3:2:n-2
        r += y[i]
    end
    r *= 2.0
    for i in (1,n)
        r += y[i]
    end
    r * dx / 3.0 #TODO: this forces the output to be Float64
end

function simps(y, x)
    n = length(y)
    if n != length(x)
        throw(ArgumentError("Input x,y must be of same length"))
    end
    if iseven(n)
        throw(ArgumentError("Simpson rule requires ODD length input (EVEN number of intervals)"))
    end

    # this is a quick generalization of the simpson's rule for
    # arbitrary separations in x
    r = zero(zero(eltype(x))*zero(eltype(y)))
    for i in 3:2:n
        d1 = x[i-1] - x[i-2]
        d2 = x[i] - x[i-1]
        h = x[i] - x[i-2]
        r += h * (
                y[i-2] * (2 - d2 / d1) +
                y[i-1] * h * h / (d1 * d2) +
                y[i] * (2 - d1 / d2)
            )
    end
    r / 6.0 #TODO: this forces the output to be Float64
end
