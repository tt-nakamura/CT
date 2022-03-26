import numpy as np
import matplotlib.pyplot as plt
from FastCT import scan_exp, reconstruct_exp

e = 4
A = plt.imread('fig1.png')[:,:,0]
A = scan_exp(A,e)
A = reconstruct_exp(A,e)
plt.imsave('fig6.png', A, cmap='gray')
