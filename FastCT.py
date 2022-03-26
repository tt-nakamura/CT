# Dicrete Radon Transform
# referece:
#   M. L. Brady, "A Fast Discrete Approximation Algorithm for the
#     Radon Transform" SIAM Journal on Computing 27 (1998) 107
#   W. H. Press, "Dicrete Radon Transform has an Exact, Fast Inverse..."
#     Proceedings of the National Academy of Sciences 103 (2006) 19249

import numpy as np
from scipy.interpolate import interpn
from CT import filtering

def scan_(A):
    """ recursive Radon transform
    input: A = image data (shape(2n,m,...))
    return B[i,j] = sum_{y=0}^{m-1} A[x,y]
             for 0<=i<2n and 0<=j<m where
           (x,y) moves from (i,0) to (i-j,m-1)
    """
    if A.shape[1] == 1: return A
    h = A.shape[1]>>1
    L,R = scan_(A[:,:h]), scan_(A[:,h:])

    A = np.repeat(L,2,1)
    for j in range(A.shape[1]):
        k = j+1>>1
        A[k:,j] += R[:A.shape[0]-k, j>>1]

    return A

def scan(A):
    """ fast Radon transform
    input: A = image data (shape(n,n,...))
    return B(r,theta) (shape(2n,n,4,...))
      B[i,j,0] = sum_{y=0}^{n-1} A[x,y] (0<=theta<=45)
        (x,y) moves from (i,0) to (i-j,n-1)
      B[i,j,1] = sum_{x=0}^{n-1} A[x,y] (45<=theta<=90)
        (x,y) moves from (0,i) to (n-1,i-j)
      B[i,j,2] = sum_{x=n-1}^0 A[x,y] (90<=theta<=135)
        (x,y) moves from (n-1,i) to (0,i-j)
      B[i,j,3] = sum_{y=n-1}^0 A[x,y] (135<=theta<=180)
        (x,y) moves from (i,n-1) to (i-j,0)
    """
    n = A.shape[0]
    if A.shape[1] != n:
        raise RuntimeError('image must be square')
    if n&(n-1):
        raise RuntimeError('n must be power of 2')

    B = np.swapaxes(A,0,1)
    C = np.zeros((n*2, n, 4) + A.shape[2:])
    C[:n,:,0] = A
    C[:n,:,1] = B
    C[:n,:,2] = B[:,::-1]
    C[:n,:,3] = A[::-1]

    return scan_(C)
   
def BackScan(A):
    """ inverse fast Radon transform
    input: A = scanned and filtered data (shape(2n,n,4,...))
    return image B restored from A (shape(n,n,...))
      B[i,j] = sum_{y=0}^{n-1} (
           A[x,y,0] (x,y) moves from (i,0) to (i+j,n-1)
         + A[x,y,1] (x,y) moves from (j,0) to (j+i,n-1)
         + A[x,y,2] (x,y) moves from (j,0) to (j+i',n-1)
         + A[x,y,3] (x,y) moves from (i',0) to (i'+j,n-1)
      )/4/(n-1)    where i'=n-1-i
    """
    n = A.shape[1]
    B = A.copy()
    B[:,[0,-1],1::2] = 0 # avoid double couting rays
    B = scan_(B[::-1])[::-1][:n]
    C = np.swapaxes(B,0,1)
    return (B[:,:,0] + B[::-1,:,3]
          + C[:,:,1] + C[::-1,:,2])/4/(n-1)

def reconstruct(A):
    return BackScan(filtering(A))

def stitch(A):
    """
      /|\   /|\
     / | \ / | \
    |  |  |  |  |
    |  |  |  |  |
    | / \ | / \ |
    |/   \|/   \|
    0 45 90 135 180
    """
    n = A.shape[1]
    n2,n4 = n<<1,n<<2
    B = np.zeros((n2,n4) + A.shape[3:])
    for j in range(n):
        k = (n-j)>>1
        B[k:,[j,n2-1-j,n2+j,-1-j]] = A[:n2-k,j]
    return B

def RadonFromSinogram(A, n=None):
    """
    input: A = output of scan() in CT.py
           n = image size
           if n is None, n is set to A.shpae[0]/2
    return: B = input to BackScan() in FastCT.py
    """
    N,M = A.shape[:2]
    if n is None: n = N>>1
    R = (n-1)/np.sqrt(2)
    r_ = np.linspace(-R,R,N)
    th_ = np.linspace(0, np.pi, M, endpoint=False)
    i,j = np.arange(n*2), np.arange(n)
    th = np.arctan2(j,n-1)
    cth = np.cos(th)
    r = (i[:,np.newaxis] - (j+n-1)/2)*cth
    th = np.c_[th, np.pi/2-th, np.pi/2+th, np.pi-th]
    p = np.broadcast_arrays(r[...,np.newaxis], th)
    p = np.moveaxis(p,0,-1)
    B = interpn((r_,th_), A, p, bounds_error=False, fill_value=0)
    B[:,0,3] = B[::-1,0,0]
    return np.einsum('ij...,j->ij...', B, cth)

def SinogramFromRadon(A, N=None, M=None):
    """
    input: A = output of scan() in FastCT.py
           N = number of parallel X-rays
           M = number of directions of X-rays
           if N is None, N is set to A.shape[1]*2
           if M is None, M is set to A.shape[1]*4
    output: B = input to BackScan() in CT.py
    """
    n = A.shape[1]
    if N is None: N = n<<1
    if M is None: M = n<<2
    R = (n-1)/np.sqrt(2)
    i,j = np.arange(n*2),np.arange(n)
    r = np.linspace(-R,R,N)
    th = np.linspace(0, np.pi, M, endpoint=False)
    B = np.empty((N,M) + A.shape[3:])
    q = np.arange(5)*np.pi/4
    for k in range(4):
        l = np.logical_and(th >= q[k], th < q[k+1])
        th1 = (th[l] - np.pi/2*(k+1>>1))*(1 - 2*(k&1))
        sc = 1/np.cos(th1)
        y = (n-1)*np.tan(th1)
        x = r[:,np.newaxis]*sc + (y+n-1)/2
        p = np.broadcast_arrays(x,y)
        p = np.moveaxis(p,0,-1)
        z = interpn((i,j), A[:,:,k], p,
                    bounds_error=False, fill_value=0)
        B[:,l] = np.einsum('ij...,j->ij...', z, sc)

    return B

############################################################
# codes below are experimental

def scan_exp(A,e):
    """ improve image quality by expansion and contraction
    e = expansion factor (power of 2)
    reference: M. L. Brady, section 3.4
    """
    n = A.shape[0]
    if A.shape[1] != n:
        raise RuntimeError('image must be square')
    if e&(e-1):
        raise RuntimeError('e must be power of 2')

    m = n*e
    B = np.zeros((m,m) + A.shape[2:])
    B[:,::e] = np.repeat(A,e,0)
    A = scan(B)
    j = ((m-1)*np.arange(n)*2 + n-1)//(2*(n-1))
    return A[::e,j]

def BackScan_exp(A,e):
    """ improve image quality by expansion and contraction
    e = expansion factor (power of 2)
    reference: M. L. Brady, section 3.4
    """
    if e&(e-1):
        raise RuntimeError('e must be power of 2')

    n = A.shape[1]
    m = n*e
    B = np.zeros((m*2,m,4) + A.shape[3:])
    B[:,::e] = np.repeat(A,e,0)
    A = BackScan(B)
    j = ((m-1)*np.arange(n)*2 + n-1)//(2*(n-1))
    return A[::e,j]

def reconstruct_exp(A, e):
    return BackScan_exp(filtering(A), e)
