import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("../main/medie_sigma_corr.txt")

x = data[:, 0]
y = data[:, 1]
erry = data[:, 2]

plt.figure()
plt.plot(x, y, erry,'o')
plt.plot(x, erry)
plt.xlabel("")
plt.ylabel("")
plt.grid(True)
plt.savefig("medie_errori_corr.png", dpi=150)
plt.show()
