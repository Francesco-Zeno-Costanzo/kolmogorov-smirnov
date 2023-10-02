#=
Code to test the quadgk routine.
We only use 4 functions but we can extend the list if we want.
The comparison is made by calculating the integrals with wolfram
=#

include("quadGK.jl")

function h(x)
    return sin(x .^2)
end
    
function f(x, a1)
    return exp(-0.5 .* (x .- a1) .^2)./sqrt(2*pi)
end
    
function g(x, a1, a2)
    return exp(-(x - a1)/a2)/(1 + exp(-(x - a1)/a2))
end

function l(x)
    num = sign(x) * atan(x^2) - abs(x)
    den = sin(x)^2 + abs(cos(x))
    return num/den
end

#=========================================================#

# should be \sim 0.7726517126900
I, dI = quadgk(h, 0, pi)
println("Integral value is ", I, " +- ", dI)

# should be 1
I, dI = quadgk(f, -100, 100, args=(1,))
println("Integral value is ", I, " +- ", dI)

# should be 10
I, dI = quadgk(g, 0, 20, args=(10, 0.02))
println("Integral value is ", I, " +- ", dI)

# should be \sim --8.734967877389
I, dI = quadgk(l, -pi, pi)
println("Integral value is ", I, " +- ", dI)
