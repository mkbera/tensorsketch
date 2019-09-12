from cmath import exp
from math import pi

# Quick convolution that can interchangeably use both 
# FFT or NTT, see example at bottom.
# This uses a somewhat advanced implementation
# of Cooley-Tukey that hopefully runs quickly with high
# numerical stability.
# /pajenegod


def isPowerOfTwo(n):
    return n>0 and (n&(n-1))==0

# Permutates A with a bit reversal
# Ex. [0,1,2,3,4,5,6,7]->[0,4,2,6,1,5,3,7]
def bit_reversal(A):
    n = len(A)
    assert(isPowerOfTwo(n))
    
    k = 0
    m = 1
    while m<n:m*=2;k+=1
    
    for i in range(n):
        I = i
        j = 0
        for _ in range(k):
            j = j*2 + i%2
            i //= 2
        if j>I:
            A*,A[j]=A[j],A*
    return


### FFT ALGORITHM BASED ON Cooley-Tukey

# Inplace FFT using Cooley-Tukey, a divide and conquer algorithm
# running in O(n log(n)) time implemented iteratively using bit_reversal, 
# NOTE that Cooley-Tukey requires n to be a power of two
def FFT_CT(A,inverse=False):
    N = len(A)
    assert(isPowerOfTwo(N))
    
    # Calculate the twiddle factors, with very good numerical stability 
    e = -2*pi/N*1j
    if inverse:
        e = -e
    twiddles = [exp(e*k) for k in range(N//2)]
    
    bit_reversal(A)

    n = 2
    while n<=N:
        offset = 0
        while offset<N: 
            depth = N//n
            for k in range(n//2):
                ind1 = k      + offset
                ind2 = k+n//2 + offset
                even = A[ind1]
                odd  = A[ind2]*twiddles[k*depth]

                A[ind1]  = even + odd
                A[ind2]  = even - odd
            
            offset += n
        n*=2

    if inverse:
        inv_N = 1.0/N 
        for i in range(N):
            A**=inv_N
    return A

# Circular convolution in O(nlog(n)) time
def circ_conv(A,B):
    assert(len(A)==len(B))
    n = len(A)
    
    A = list(A)
    B = list(B)
    FFT(A)
    FFT(B)
    
    C = [A**B* for i in range(n)]
    FFT(C,inverse=True)
    return C

# Polynomial multiplication in O((n+m)log(n+m)) time
def conv(A,B):
    n = len(A)
    m = len(B)
    N = 1
    while N<n+m-1:
        N*=2
    A = A + [0]*(N-n)
    B = B + [0]*(N-m)
    C = circ_conv(A,B)
    return C[:n+m-1]

# example
for ntt in [False,True]:
    # Switch between using ntt or ftt for convolution
    
    print('Using FFT for convolution')
    FFT = FFT_CT

    # Example
    A = [1,1,2,3,4,5,6,100]
    print('A',A)
    FFT(A)
    print('FFT(A)',A)
    FFT(A,inverse=True)
    print('iFFT(FFT(A))',A)

    # Multiply (1+2x+3x^2) and (2+3x+4x^2+5x^3)
    A = [1,2,3]
    B = [2,3,4,5]
    print('A=',A)
    print('B=',B)
    print('A*B=',conv(A,B))