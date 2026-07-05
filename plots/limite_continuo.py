import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def fit_polinomiale_3deg(x, c0, c1, c2, c3):
    return c0 + c1*x + c2*x**2 + c3*x**3

def f_true(x):
    return 1 - (1/24)*x + (3/640)*x**2 - (5/7168)*x**3

# Carica i dati
data = np.loadtxt("dati.txt")

a    = data[:, 0]        # passo reticolare
x_a2 = a**2             # calcola a²
y    = data[:, 1]        # Delta E
erry = data[:, 2]        # errore

popt, pcov = curve_fit(
    fit_polinomiale_3deg,
    x_a2,
    y,
    sigma=erry,
    absolute_sigma=True
)

c0, c1, c2, c3 = popt
perr = np.sqrt(np.diag(pcov))

residuals = y - fit_polinomiale_3deg(x_a2, *popt)
chi2      = np.sum((residuals / erry)**2)
ndof      = len(x_a2) - len(popt)
chi2_red  = chi2 / ndof

print("=========================================")
print("RISULTATI FIT:")
print("=========================================")
print(f"c0 = {c0:.6f} ± {perr[0]:.6f}  (atteso: 1.0)")
print(f"c1 = {c1:.6f} ± {perr[1]:.6f}  (atteso: -0.04167)")
print(f"c2 = {c2:.6f} ± {perr[2]:.6f}  (atteso:  0.00469)")
print(f"c3 = {c3:.6f} ± {perr[3]:.6f}  (atteso: -0.00070)")
print(f"\nchi2 ridotto = {chi2_red:.2f}  (ndof = {ndof})")
print("=========================================")

# Grafico
x_plot = np.logspace(np.log10(min(x_a2)*0.5), np.log10(max(x_a2)*1.5), 500)

plt.figure(figsize=(8, 5))
plt.errorbar(x_a2, y, yerr=erry, fmt='o', color='tab:blue',
             capsize=4, label='Dati')
plt.plot(x_plot, fit_polinomiale_3deg(x_plot, *popt), '-',
         color='tab:blue', label='Fit')
plt.plot(x_plot, f_true(x_plot), '--',
         color='tab:orange', label='Teoria')
plt.axhline(1.0, color='gray', linestyle=':', linewidth=0.8)
plt.xscale('log')
plt.xlabel(r'$a^2$')
plt.ylabel(r'$\Delta E$')
plt.title('Limite al continuo')
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig("fit_continuo.png", dpi=150)
plt.show()
