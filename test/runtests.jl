using Quadrature
using Base.Test

xsamp = linspace(0.0, 10.0, 11)
ysamp = map(x->x^3, xsamp)
@test simps(ysamp, xsamp) == 2500.0

xsamp = linspace(0.0, 10.0, 1001)
ysamp = map(x->x^4, xsamp)
@test_approx_eq simps(ysamp, xsamp) 20000.0

# make sure we give errors for wrong input shape
# input ans samples are not the same length
@test_throws ArgumentError simps([0.5, 0.5, 0.3], [0.1, 0.2])
# simpsons rule needs odd input length
@test_throws ArgumentError simps([0.5, 0.5, 0.3, 0.2], [0.1, 0.2, 0.3, 0.4])
