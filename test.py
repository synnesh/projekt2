import numpy as np
import sys


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


def jakobi_rotate(A,R,k,l,n):

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

	A[k,k] = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s
	A[l,l] = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s

	A[k,l] = 0
	A[l,k] = 0

	for i in range(n):
		if (i != k and i != l):
			a_ik = A[i,k]
			a_il = A[i,l]
			A[i,k] = a_ik*c - a_il*s
			A[i,l] = a_il*c + a_ik*s
			A[k,i] = A[i,k]
			A[l,i] = A[i,l]

        r_ik = R[i,k];
        r_il = R[i,l];

        R[i,k] = c*r_ik - s*r_il;
        R[i,l] = c*r_il + s*r_ik;
	return A, R


def test_offdiag():
    A = np.matrix([[1,7,1,1],[7,1,1,1],[1,1,1,1],[1,1,1,1]])
    expected_p = 0
    expected_q = 1
    computed_p,computed_q = offdiag(A,4)
    tol = 1E-10
    success = abs(computed_p-expected_p < tol and computed_q-expected_q < tol)
    assert success

def test_offdiagmaks():
    A = np.matrix([[1,7,1,1],[7,1,1,1],[1,1,1,1],[1,1,1,1]])
    expected = 7
    computed = offdiagmaks(A,4)
    tol = 1E-10
    success = abs(computed - expected < tol)
    assert success


toleranse = 1.0e-10
N = 0
maxiter = 10000
maks = 1;
n= 50
ro_min = 0
ro_max = 5
R = np.identity(n)
testiter=1


def matrise(ro_min, ro_max, n):
    h = float(ro_max - ro_min)/(n+1)
    A = np.zeros((n,n))

    e = -1./(h**2)
    for i in range(n):
        Vi = (ro_min + i*h)**2
        A[i,i]=2./(h**2)+Vi

    for j in range(n-1):
        A[j,j+1] = e
        A[j+1,j] = e
    return A

A = matrise(ro_min,ro_max,n)
print np.sort(np.linalg.eigvals(A))

while(offdiagmaks(A,n) > toleranse and N <= maxiter):
    k,l = offdiag(A,n)
    A,R = jakobi_rotate(A,R,k,l,n)
    test = N/100.
    tol = 10E-10
    if (abs(test-testiter)<tol):
        testiter+=1
        if (np.inner(R[1,:],R[5,:]) > tol):
            print "The orthogonality is not preserved for the eigenvectors"
            sys.exit([1])
    N+=1
print testiter
print np.sort(A.diagonal())

print N
test_offdiag()
test_offdiagmaks()
