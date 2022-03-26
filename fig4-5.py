import numpy as np
import matplotlib.pyplot as plt
from FastCT import scan,stitch,reconstruct

A = plt.imread('fig1.png')[:,:,0]
A = scan(A)

# fig4
r = np.r_[0:1:5j]*A.shape[0]
th = np.r_[0:1:5j]*A.shape[0]*2
plt.figure(figsize=(6,3.2))
plt.imshow(stitch(A), cmap='gray')
plt.xticks(th, ('0','45','90','135','180'))
plt.yticks(r, ('-2','-1','0','1','2'))
plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
plt.close()

# fig5
A = reconstruct(A)
plt.imsave('fig5.png', A, cmap='gray')
