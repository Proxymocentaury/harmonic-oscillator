import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("../main/termalizzazione.txt")

sweep  = data[:, 0]
action = data[:, 1]

plt.figure()
plt.plot(sweep, action)
plt.xlabel("sweep")
plt.ylabel("S")
plt.title("Termalizzazione")
plt.grid(True)
plt.savefig("termalizzazione.png", dpi=150)
plt.show()
