import numpy as np
import matplotlib.pyplot as plt 

t = np.linspace(0,50,1000)
g = np.exp(-(t-3-1)/2) - np.exp(-(t-3-1)/0.4)

plt.plot(t,g, 'k.-')
plt.xlim(4, 50)
plt.ylim(0, 1)
plt.show()