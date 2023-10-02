#=
Code to test the ran2 routine.
=#

# Must be included in the REPL to speed up
# using Plots
include("ran2.jl")
ranstart("randomseed")

function gauss(x, mu, s)
    # normal distribution
    g = exp.(-0.5 .* (x .- mu) .^2 /s .^2)./sqrt(2*pi*s^2)
    return g
end

#==============================
Uniform test
==============================#

N = Int(1e6)
y = [ran2() for i in 1:N]

println("computed vs analytical")
for i in 1:10
    m_s = sum(y .^i)/N # extimeted momentum
    m_t = 1/(1+i)      # analytical momentum
    println(i, "-th momentum: ", m_s, " ",  m_t)
end

p1 = histogram(y, bins=50, normalize=:pdf, label="uniform", minorgrid=true)
title!("ran2 distribution")
xlabel!("x")
ylabel!("p(x)")

#==============================
Normal test
==============================#

x1 = ran2(N)
x2 = ran2(N)

s = 1
m = 0

# BOX MULLER 
phi = 2 .* pi .* x1
z = -log.(1 .- x2)
y1 = sqrt.(2 .* z .* s^2) .*cos.(phi) .+ m
y2 = sqrt.(2 .* z .* s^2) .*sin.(phi) .+ m

p2 = histogram(y1, bins=50, normalize=:pdf, label="Normal", yticks=0:0.05:0.4)
t = -5:0.1:5
p3 = plot(p2, t, gauss(t, m, s), lc=:black, label="gaussian", minorgrid=true)
title!("box muller")
xlabel!("x")
ylabel!("p(x)")

ranfinish("randomseed")

plot(p1, p3)
