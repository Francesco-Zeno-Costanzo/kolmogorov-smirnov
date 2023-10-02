
int(x) = floor(Int32, x) # int function, e.g. int(2.3) = 2

#=============================================================
Parameter for ran2 generator
=============================================================#
IM1  = 2147483563
IM2  = 2147483399
AM   = 1.0/IM1
IMM1 = IM1 - 1
IA1  = 40014
IA2  = 40692
IQ1  = 53668
IQ2  = 52774
IR1  = 12211
IR2  = 3791
NTAB = 32
NDIV = 1 + int(IMM1/NTAB)
EPS  = 1.2e-7
RNMX = 1.0 - EPS # largest floating value that is less than 1
#=============================================================
Parameter for ranstart and ranfinish for update seed file
=============================================================#
idum  = 0
idum2 = 0
iv    = [0 for l in 1:32]
iy    = 0
#=============================================================
RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
=============================================================#


@doc raw"""
Function to reed seed file
    
Parameter
---------
path : string
    path of seed file
"""
function ranstart(path)

    # global variables
    global idum, idum2, iv, iy, NTAB

    f = open(path)
    # read the first value, i.e. the initial seed
    idum  = parse(Int32, readline(f))
    
    # in case ranfinish has already been called previously
    try   # try to reed the second one
        idum2 = parse(Int32, readline(f))
    catch # if nothing is written I set idum=0
        idum2 = 0
    else  # else use read value
        idum2 = idum2
    end
    
    if idum2 > 0 
        for i in 1:NTAB
            iv[i] = parse(Int32, readline(f))
        end
        iy = parse(Int32, readline(f))
    else
        # in the case we are at the first call
        if idum >= 0
            idum = - idum - 1
        end
    end
    
    close(f)
    
end


@doc raw"""
Returns a uniform random variable within (0.0, 1.0)
The idea of ran2 is to combine two different sequences with
different periods so as to obtain a new sequence whose period
is the least common multiple of the two periods. We use :

    I_{i+1} = a*I_i mod m; with m = a*q + r 
    The parameters of the two generators are:
    m1 = 2147483563, a1 = 40014, q1 = 53668, r1 = 12211
    m2 = 2147483399, a2 = 40692, q2 = 52774, r2 = 3791

Both generators' moduli are slightly less than 2^31; and their
factorization is: 

    m1 − 1 = 2 × 3  × 7  × 631  × 81031 
    m2 − 1 = 2 × 19 × 31 × 1019 × 1789 

so the period of ran2 is \simeq 2.3 x 10^{18}

Parameters
----------
size : Nothing, optional, defaul nothing \
    if you call ran2() you pass nothing and
    only one number will be generated
        
Return
------
pseudo-random number

Usage
-----
To be used you need to create a file called with a number
that will act as a seed. The file will be read by the
ranstart function and updated by the ranfinish function.
Writing what will be necessary for subsequent generations.
e.g.:

# Start of the code:
    path = ... # path of seed file
    ranstart(path) # reed the seed

#Calculations in between: 
    x = ran2()     # generate a random number

# It is also possible to generate an array of length N thanks to multiple dispatch 
    N = ...        # size of array
    y = ran2(N)    # generate an array of random number 

# End of the code
    ranfinsh(path) # update seed for future use

"""
function ran2(size::Nothing=nothing)

    # global variables
    global IM1, IM2, AM, IMM1, IA1, IA2
    global IQ1, IQ2, IR1, IR2, NTAB, NDIV 
    global EPS, RNMX, idum, idum2, iv, iy 
    
    if idum <= 0                          # initialize
        idum = max(-idum, 1)              # prevent idum = 0
        idum2 = idum
        for j in NTAB+8:-1:1              # load shuflle table
            k = int(idum/IQ1)
            idum = IA1*(idum - k*IQ1) - k*IR1
            if idum < 0
                idum = idum + IM1
            end
            if j <= NTAB
                iv[j] = idum
            end
        end
        iy = iv[1]
    end
    k = int(idum/IQ1)                     # Start here when not initializing.
    idum = IA1*(idum - k*IQ1) - k*IR1     # Compute idum=mod(IA1*idum,IM1) 
    if idum < 0                           # without overflows by Schrage’s method.
        idum = idum + IM1
    end
    k = int(idum2/IQ2)
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2   # Compute idum2=mod(IA2*idum2,IM2)
    if idum2 < 0
        idum2 = idum2 + IM2
    end
    j = 1 + int(iy/NDIV)                  # Will be in the range 1:NTAB.
    iy = iv[j] - idum2                    # Here idum is shuffled, idum and idum2
    iv[j] = idum                          # are combined to generate output.
    if iy < 1
        iy = iy + IMM1
    end
    return min(AM*iy, RNMX) #Because users don’t expect endpoint values
    
end


@doc raw"""
Returns an array of random number within (0.0, 1.0)
according a uniform distribution.
    
Parameters
----------
size : Int
    size of the array that must be generated.
    The generation is done by calling N time ran2(). 
    
Return
------
R : array
    Returns an array of random number within (0.0, 1.0)    
"""
function ran2(size::Int)    
    R = [ran2() for i in 1:size]
    return R
end


@doc raw"""
Function to update seed file
    
Parameter
---------
path : string
    path of seed file
"""
function ranfinish(path)
    
    # global variables
    global idum, idum2, iv, iy, NTAB
    
    f = open(path, "w")
    write(f, string(idum))
    write(f, "\n")
    write(f, string(idum2))
    write(f, "\n")
    for i in 1:NTAB
        write(f, string(iv[i]))
        write(f, "\n")
    end
    write(f, string(iy))
    write(f, "\n")
    close(f)
end
