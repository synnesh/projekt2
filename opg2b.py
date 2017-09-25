import numpy as np


def offdiag(A, n):
    maks = 0
    for i in range(n):
        for j in range(i+1,n):
          aij = abs(A[i,j])
          #print aij
          if (aij > maks):
              maks = aij
              p = i
              q = j
    return p, q

def offdiagmaks(A, n):
    maks = 0
    for i in range(n):
        for j in range(i+1,n):
          aij = abs(A[i,j])
          #print aij
          if (aij > maks):
              maks = aij
    return maks


def jakobi_rotate(A,k,l,n):

    if(A[k,l] != 0):
        tau = float(A[l,l] - A[k,k])/(2*A[k,l])

        if (tau >= 0):
            t = 1.0/(tau + np.sqrt(1.0 + tau*tau))
        else:
            t = -1.0/(-tau + np.sqrt(1.0 + tau*tau))
        c = 1.0/np.sqrt(1+t*t)
        s = c*t
    else:
        c = 1.0
        s = 0.0


    a_kk = A[k,k]
    a_ll = A[l,l]
    a_kl = A[k,l]

    A[k,k] = a_kk*c**2 - 2*a_kl*c*s + a_ll*s**2
    A[l,l] = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s

    A[k,l] = 0
    A[l,k] = 0

    for i in range(n):
        if (i != k and i != l):
            A[i,k] = A[i,k]*c - A[i,l]*s
            A[i,l] = A[i,l]*c + A[i,k]*s
            A[k,i] = A[i,k]
            A[l,i] = A[i,l]

    return A

toleranse = 1.0e-10
N = 0
maxiter = 1000
maks = 1;
n= 40
ro_min = 0
ro_max = 10

def matrise(ro_min, ro_max, n):
    h = float(ro_max - ro_min)/n
    A = np.zeros((n,n))
    e = -1./h**2
    for i in range(n):
        Vi = (ro_min + i*h)**2
        A[i,i]=2./(h**2)+Vi
    for j in range(n-1):
        A[j,j+1] = e
        A[j+1,j] = e
    return A

A = matrise(ro_min,ro_max,n)
print np.linalg.eig(A)
#print A
while(offdiagmaks(A,n) > toleranse and N <= maxiter):
    k,l = offdiag(A,n)
    A = jakobi_rotate(A,k,l,n)
    N+=1
print A.diagonal()
