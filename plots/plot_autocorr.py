import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --------------------
# Lettura dati
# --------------------

autocorr     = np.loadtxt("../main/autocorrelatore.txt")
autocorr_bin = np.loadtxt("../main/autocorrelatorebin.txt")

# --------------------
# TAGLIO INIZIALE (qui decidi il range)
# --------------------

t_max = 100

autocorr     = autocorr[autocorr[:, 0] <= t_max]
autocorr_bin = autocorr_bin[autocorr_bin[:, 0] <= t_max]

# separazione dopo il taglio
t1 = autocorr[:, 0]
g1 = autocorr[:, 1]

t2 = autocorr_bin[:, 0]
g2 = autocorr_bin[:, 1]

# --------------------
# Fit esponenziale (solo primo dataset)
# --------------------

def expo(t, A, tau):
    return A * np.exp(-t / tau)

# fit su dati già filtrati
mask = g1 > 0

popt, pcov = curve_fit(
    expo,
    t1[mask],
    g1[mask],
    p0=(g1[0], 5.0)
)

A_fit, tau_fit = popt

tfit = np.linspace(0, t_max, 300)
gfit = expo(tfit, A_fit, tau_fit)

# --------------------
# Plot
# --------------------

plt.figure(figsize=(8,6))

plt.plot(t1, g1, 'o', label="Autocorrelatore")
plt.plot(t2, g2, 'o', label="Autocorrelatore binnato")

plt.plot(tfit, gfit, '--', label=rf"Fit exp, $\tau={tau_fit:.3f}$")

plt.axhline(0, color='black', linestyle='--', linewidth=0.8)

plt.xlim(0, t_max)

plt.xlabel("t")
plt.ylabel(r"$\Gamma(t)$")
plt.title("Autocorrelatore")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("autocorrelatore.png", dpi=150)
plt.show()