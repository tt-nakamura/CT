import numpy as np
import matplotlib.pyplot as plt
import CT
from FastCT import RadonFromSinogram, reconstruct

try:
    A = np.load('fig2.npy')
except IOError:
    A = plt.imread('fig1.png')[:,:,0]
    A = CT.scan(A)
    np.save('fig2.npy', A)

A = RadonFromSinogram(A)
A = reconstruct(A)
plt.imsave('fig7.png', A, cmap='gray')
