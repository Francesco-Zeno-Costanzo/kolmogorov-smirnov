#=
test of KS tets
=#
#using PythonPlot # include first in repl
using Distributions
include("KS_dist.jl")
include("ran2.jl")
ranstart("randomseed")


function gauss(x, mu, s)
    # gaussian distribution
    g = exp.(-0.5 .* (x .- mu) .^2 /s .^2)./sqrt(2*pi*s^2)
    return g
end

#==============================================================
Sample generation
==============================================================#
M = 3000
s0 = 1
m0 = 0
s1 = 0.8
m1 = 3
# uniform
x1 = ran2(M)
x2 = ran2(M)
# BOX MULLER 
phi = 2 .* pi .* x1
z = -log.(1 .- x2)
y1 = sqrt.(2 .* z .* s0^2) .*cos.(phi) .+ m0
y2 = sqrt.(2 .* z .* s1^2) .*sin.(phi) .+ m1

#==============================================================
Test distribution: P value of set of measurament
==============================================================#

KS1, p_value1 = sf(y1, gauss, -15.; args=(m0, s0), N=Int(1e4))
println(KS1, " ", p_value1) # accept H_0
KS2, p_value2 = sf(y2, gauss, -15.; args=(m0, s0), N=Int(1e4))
println(KS2, " ", p_value2) # reject H_0

t = -5:0.01:8
fig1 = figure(1)
hist(y1, bins=50, density=true, label="data1", histtype="step")
hist(y2, bins=50, density=true, label="data2", histtype="step")
plot(t, gauss(t, m0, s0), 'k', label="N($(m0), $(s0))")
plot(t, gauss(t, m1, s1), 'b',  label="N($(m1), $(s1))")
title("Data distributions")
xlabel("x")
ylabel("p(x)")
grid()

#==============================================================
Test distribution: cfr with distribution library
==============================================================#

# from code
KS = KS_compute_dist(M)
XKS, KSCDF = compute_cdf(KS)

# from Distributions
N     = 1000
x_d   = maximum(KS)
dx    = x_d/(N-1)
xks   = Array(0:dx:x_d)
ks    = KSDist(M)
kscdf = cdf.(ks, xks)
# compute pdf as derivative of cdf, because pdf is not implemented
dcdf = zeros(N)
dcdf[1]     = (- 3*kscdf[1]  + 4*kscdf[2]  - kscdf[3] )/(2*dx)
dcdf[2:N-1] = (kscdf[3:N] - kscdf[1:N-2])/(2*dx)
dcdf[N]     = ( 3*kscdf[N] - 4*kscdf[N-1] + kscdf[N-2])/(2*dx)

fig2 = figure(2)
subplot(121)
hist(KS, bins=50, density=true, label="KS", histtype="step")
plot(xks, dcdf, 'k', label="Large N")
title("KS distribution")
xlabel("x")
ylabel("p(x)")
grid()

subplot(122)
plot(xks, 1 .- kscdf, 'b', label="Large N")
plot(XKS, 1 .- KSCDF, 'k', label="KS")
title("Survival function")
xlabel("x")
ylabel("sv(x)")
grid()

ranfinish("randomseed")
