import numpy as np
import matplotlib.pyplot as plt
from CT import scan,BackScan,filtering

try:
    A = np.load('fig2.npy')
except IOError:
    A = plt.imread('fig1.png')[:,:,0]
    A = scan(A)
    np.save('fig2.npy', A)

# fig2
a = np.sqrt(2)/2
r = np.r_[1-a,1,1+a]*A.shape[0]/2
th = np.r_[0:1:5j]*A.shape[1]
plt.figure(figsize=(5.5, 6))
plt.subplot(2,1,1)
plt.imshow(A, cmap='gray')
plt.xticks(th, ('0','45','90','135','180'))
plt.yticks(r, ('-1','0','1'))
plt.ylabel(r'$r$')
plt.subplot(2,1,2)
A = filtering(A)
plt.imshow(A, cmap='Greys')
plt.xticks(th, ('0','45','90','135','180'))
plt.yticks(r, ('-1','0','1'))
plt.xlabel(r'$\theta$')
plt.ylabel(r'$r$')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
plt.close()

# fig3
A = BackScan(A)
plt.imsave('fig3.png', A, cmap='gray')
