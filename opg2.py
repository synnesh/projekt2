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
    return p, q, maks

def offdiagmaks(A, n):
    maks = 0
    for i in range(n):
        for j in range(i+1,n):
          aij = abs(A[i,j])
          #print aij
          if (aij > maks):
              maks = aij
    return maks



def jakobi_rotate(A, R, k, l, n):
    if (A[k,l] != 0):
        tau = (A[l,l] - A[k,k])/(2*A[k,l])

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
    A[k,k] = c*c*a_kk - 2.0*c*s*A[k,l] + s*s*a_ll;
    A[l,l] = s*s*a_kk + 2.0*c*s*A[k,l] + c*c*a_ll;
    A[k,l] = 0.0;  # hard-coding non-diagonal elements by hand
    A[l,k] = 0.0;  # same here
    for i in range(n):
        if ( i != k and i != l ):
            a_ik = A[i,k]
            a_il = A[i,l]
            A[i,k] = c*a_ik - s*a_il
            A[k,i] = A[i,k]
            A[i,l] = c*a_il + s*a_ik
            A[l,i] = A[i,l]
##  And finally the new eigenvectors
        r_ik = R[i,k];
        r_il = R[i,l];

        R[i,k] = c*r_ik - s*r_il;
        R[i,l] = c*r_il + s*r_ik;
        return A,R



def nondiagmax(A):
    return np.fill_diagonal(A, 0).max()


A = np.matrix('6 -2 -1; -2 6 -1; -1 -1 5')
R = np.zeros((3,3))
np.fill_diagonal(R,1)
toleranse = 1.0e-10
N = 0
maxiter = 100
maks = 1;

while(offdiagmaks(A,3) > toleranse and N <= maxiter):
    a = offdiag(A,3)
    k = a[0]
    l = a[1]
    maks = a[2]
    A,R = jakobi_rotate(A,R,k,l,3)
    N+= 1
    print k,l, maks

print A, R
