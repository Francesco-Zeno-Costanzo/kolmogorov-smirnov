#=
Code for creating kolmogorov smirnov distribution.
=#
include("ran2.jl")
include("quadGK.jl")
ranstart("randomseed")

@doc raw"""
Function for compute empirical cdf
    
Parameters
----------
data : array
    array of sampled data
    
Retunrs
-------
(data_s, cdf) : tuple
    tuple with sorted data and cdf
"""
function compute_cdf(data::Vector)
    n = length(data)
    data_s = sort(data)
    cdf = [i/n for i in 1:n]
    return (data_s, cdf)
end


@doc raw"""
Function for compute numerical cdf
    
Parameters
----------
func : function
    pdf to integrate \
range : tuple
    range of intergral \
args : tuple, optional
    extra agruments to pass at func
    
Retunrs
-------
I : flota
    value of integral
"""
function compute_cdf(func::Function, range::Tuple; args::Tuple=())
    I, dI = quadgk(func, range[1], range[2], args=args, tol=1e-7)
    return I  
end


@doc raw"""
Function for compute KS statistics.
We use an uniform distribution beacuse it is known 
that the distribution of KS does not depend on the
theorical pdf used and the uniform distribution
speeds up the calculation. 

Parameters
----------
data : array
    array of sampled data

Retunrs
-------
DN : float
    value of KS for input data
"""
function KS_test(data::Vector)
    data_cdf = compute_cdf(data)
    cdf = data_cdf[1] # assuming uniform distribution
    Dn = maximum(abs.(data_cdf[2] .- cdf))
    return Dn
end


@doc raw"""
Function for compute KS statistics.
Using a numerical integration for theoretical pdf
    
Parameters
----------
data : array
    array of sampled data \
f : function
    theoretical pdf \
start : float
    lower bound for integration \
args : tuple
    extra argument to pass at f
    
Retunrs
-------
DN : float
    value of KS for input data
"""
function KS_test(data::Vector, f::Function, start::Float64; args::Tuple=())
    data_cdf = compute_cdf(data)
    cdf = [compute_cdf(f, (start, x), args=args) for x in data_cdf[1]] 
    Dn = maximum(abs.(data_cdf[2] .- cdf))
    return Dn
end


@doc raw"""
Function for compute KS distribution.

Parameters
----------
M : Int
    Size of sampled data \
N : Int, optional, default Int(1e4)
    Number of sample to build the KS distribution

Retunrs
-------
KS : array of lenght N
    sampled KS  
"""
function KS_compute_dist(M::Int; N::Int=Int(1e4))
    KS = zeros(N)
    for i in 1:N
        KS[i] = KS_test(ran2(M))
    end
    return KS
end


@doc raw"""
Function for compute p-value from survival function.
If the KS value is too large and goes outside the bounds of the constructed
distribution, the p-value cannot be calculated and the smallest non-zero
value is returned, which should therefore be understood as an overestimate.

Parameters
----------
data : array
    array of sampled data \
f : function
    theoretical pdf \
start : float
    lower bound for integration \
args : tuple
    extra argument to pass at f \
N : Int, optional, default Int(1e4)
    Number of sample to build the KS distribution
    
Return
------
pv : float
    p-value \
KS : float
    Value of KS for data sample
"""
function sf(data::Vector, f::Function, start::Float64; args::Tuple=(), N::Int=Int(1e4))
    
    M   = length(data)                        # size of sample
    KS  = KS_test(data, f, start, args=args)  # value for our sample
    
    KSd = KS_compute_dist(M, N=N)             # KS distributions
    KS_sample, KS_cdf = compute_cdf(KSd)      # cdf of KS
    sv = 1 .- KS_cdf                          # survival function
    # find index to obtain p-value
    j = 0
    for i in 1:N-1
        if KS_sample[i] < KS && KS < KS_sample[i+1]
            j = i
        end 
    end
    if j == 0
        pv = sv[N-1]
    else
        pv = sv[j] 
    end
    return KS, pv
end


ranfinish("randomseed")
