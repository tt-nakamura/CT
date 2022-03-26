import numpy as np
from scipy.interpolate import interpn,interp1d
from numpy.fft import rfft,irfft,rfftfreq

def scan(A, n=None, m=None):
    """ sinogram (Radon transform) of image
    A = image data (shape(M,N,...))
    n = number of paralell X-rays (also used as
        number of integration points along an X-ray) 
    m = number of directions of X-rays
    return sinogram of A (shape(n,m,...))
    """
    M,N = A.shape[:2] # M=height, N=width
    if n is None: n = 1<<M.bit_length() # (power of 2) > M
    if m is None: m = n<<1
    X,Y = (M-1)/2, (N-1)/2
    R = np.hypot(X,Y)

    x,y = np.r_[-X:X:M*1j], np.r_[-Y:Y:N*1j]
    r,s = np.mgrid[-R:R:n*1j,-R:R:n*1j]
    # 0 <= theta < pi (direction of X-rays)
    th = np.linspace(0, np.pi, m, endpoint=False)
    B = np.empty((n,m) + A.shape[2:])
    for i,th in enumerate(th):
        p = [r*np.cos(th) - s*np.sin(th),
             r*np.sin(th) + s*np.cos(th)]
        p = np.moveaxis(p,0,-1)
        z = interpn((x,y), A, p, bounds_error=False, fill_value=0)
        B[:,i] = np.sum(z, axis=1)

    return B

def BackScan(A, M=None, N=None):
    """ inverse Radon transform of sinogram
    A = sinogram after filtering (shape(n,m,...))
    M,N = height and width of output image
    return image restored from A (shape(M,N,...))
    """
    n,m = A.shape[:2]
    if M is None: M = n>>1
    if N is None: N = M
    X,Y = (M-1)/2, (N-1)/2
    R = np.hypot(X,Y)

    x,y = np.mgrid[-X:X:M*1j,-Y:Y:N*1j]
    r = np.linspace(-R,R,n)
    th = np.linspace(0, np.pi, m, endpoint=False)
    B = np.zeros((M,N) + A.shape[2:])
    for i,th in enumerate(th):
        f = interp1d(r, A[:,i], axis=0, copy=False)
        B += f(x*np.cos(th) + y*np.sin(th))

    return B/m

def filtering(A):
    """ apply high-pass filter to sinogram """
    f = rfftfreq(A.shape[0])*np.pi
    A = rfft(A, axis=0)
    A = np.einsum('i,ij...->ij...', f, A)
    return irfft(A, axis=0)

def reconstruct(A, M=None, N=None):
    return BackScan(filtering(A), M, N)
