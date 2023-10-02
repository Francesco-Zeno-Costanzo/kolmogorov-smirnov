#=
Routine for the calculation of integrals using adaptive Gaussian quadrature.
Only finite interval are implemented
=#
include("gkdata.jl")


@doc raw"""
Calculation of the integral of f
with the Gauss-Kronrod quadrature method

Parameters
----------
f : function
    function to integrate \
a : float 
    beginning of the interval \
b : float 
    end of the interval \
args : tuple, optional 
    extra arguments to pass to f \
set : tuple, optional, default (G10, K21) 
    possible set to use, the aviable are: \
    (G7,  K15), (G10, K21), (G15, K31), \
    (G20, K41), (G25, K51), (G30, K61)

Return
------
I_K : float 
    \int_a^b f \
error : float 
    error of integration
"""
function gausskronrod(f, a, b; args=(), set=(10, 21))

    N, M = set
    # b must be grather tha a
    if b < a
        a, b = b, a
    else
        # do nothing
    end

    mid = 0.5 * (b + a)
    dx  = 0.5 * (b - a)

    I_G = 0
    I_K = 0

    for (gwi, gxi) in zip(G_weights[N], G_nodes[N])
        I_G += gwi*dx*f((gxi + 1)*dx + a, args...)
    end
    
    for (kwi, kxi) in zip(K_weights[M], K_nodes[M])
        I_K += kwi*dx*f((kxi + 1)*dx + a, args...)
    end
    
    error = abs(I_G - I_K)

    return I_K, error
end
   
 
@doc raw"""
Compute the integral of f; if tolerance is not satisfied,
we split the interval in left and right part and recompute
the integral and so on recursively until the tolerance is
satisfied.

Parameters
----------
f : function 
    function to integrate \
a : float 
    beginning of the interval \
b : float 
    end of the interval \
args : tuple, optional 
    extra arguments to pass to f \
set : tuple, optional, default (G10, K21) 
    possible set to use, the aviable are: \
    (G7,  K15), (G10, K21), (G15, K31), \
    (G20, K41), (G25, K51), (G30, K61) \
tol : float, optional 
    tollerance, default 1e-10 

Return 
------
I : float 
    \int_a^b f \
err : float 
    error on I \
"""
function quadgk(f, a, b; args=(), set=(10, 21), tol=1e-10)

    # chek for integration
    val_g = [7,  10, 15, 20, 25, 30]
    val_k = [15, 21, 31, 41, 51, 61]
    L = [string("G", g, ", K", k) for (g, k) in zip(val_g, val_k)]
    l = string("G", set[1], ", K", set[2])
    
    if l ∉ L
        error("Impossible to use $l, you must choose between:\n$(join(L, ", "))")
    end
    
    # compute the integral
    I, err = gausskronrod(f, a, b, args=args, set=set)
    # stop criteria
    if err < tol
        return I, err
    end

    # mid point
    m = a + (b - a)/2
    # recursive calls
    I1, err1 = quadgk(f, a, m, args=args, set=set, tol=tol)
    I2, err2 = quadgk(f, m, b, args=args, set=set, tol=tol)

    I = I1   + I2
    err = err1 + err2

    return I, err
end
