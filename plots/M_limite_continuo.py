import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ==========================================
# CONFIGURAZIONE
# ==========================================
nome_file = "dati.txt"  # Sostituisci con il nome del tuo file txt

# Modello teorico analitico dell'energia
def f_teorica(x_a2):
    return 0.5*(1 - (1/18)*x_a2 + (3/128)*x_a2**2 - (5/1024)*x_a2**3 )

# Modello di fit: Polinomio di 3° grado in a²
# y = c0 + c1*(a²) + c2*(a²)² + c3*(a³)³
def fit_cubico(x, c0, c1, c2):
    return c0 + c1*x + c2*x**2 

# ==========================================
# CARICAMENTO E PREPARAZIONE DATI
# ==========================================
# Carica le colonne: 0->a, 1->E, 2->Eerr
data = np.loadtxt(nome_file)

a_val = data[:, 0]
y_energia = data[:, 1]
err_energia = data[:, 2]
x_a2 = a_val**2

x_fit = x_a2
y_fit = y_energia
err_fit = err_energia

    
# ==========================================
# ESECUZIONE DEL FIT
# ==========================================
popt, pcov = curve_fit(
    fit_cubico,
    x_fit,
    y_fit,
    sigma=err_fit,
    absolute_sigma=True
)

c0, c1, c2 = popt
perr = np.sqrt(np.diag(pcov))

# Calcolo del Chi2 ridotto (CHI/dof)
residuals = y_fit - fit_cubico(x_fit, *popt)
chi2 = np.sum((residuals / err_fit) ** 2)
ndof = len(x_fit) - len(popt)
chi2_red = chi2 / ndof

# ==========================================
# STAMPA DEI RISULTATI
# ==========================================
print("\n=========================================")
print("          RISULTATI DEL FIT              ")
print("=========================================")
print(f"Intercetta al continuo (c0) = {c0:.6f} ± {perr[0]:.6f} (Atteso: 0.5)")
print(f"Coefficiente ord 1 (c1)     = {c1:.6f} ± {perr[1]:.6f} (Atteso: 1/18= 0.055555556)")
print(f"Coefficiente ord 2 (c2)     = {c2:.6f} ± {perr[2]:.6f}")
#print(f"Coefficiente ord 3 (c3)     = {c3:.6f} ± {perr[3]:.6f}")
print("-----------------------------------------")
print(f"Gradi di libertà (dof)      = {ndof}")
print(f"χ² totale                   = {chi2:.4f}")
print(f"χ² ridotto (χ²/dof)         = {chi2_red:.4f}")
print("=========================================\n")

# ==========================================
# GRAFICO
# ==========================================
# Genera punti densi per disegnare le curve fluide (in scala logaritmica)
x_plot = np.logspace(np.log2(min(x_a2)*0.5), np.log2(max(x_a2)*1.2), 500, base=2)

plt.figure(figsize=(9, 6))

# Disegna TUTTI i punti Monte Carlo con le barre d'errore
plt.errorbar(x_a2, y_energia, yerr=err_energia, fmt='o', color='tab:blue', 
             capsize=4, markersize=6, label='Dati Monte Carlo (10^7 conf)')

# Disegna la curva di Fit e la curva Teorica asintotica
plt.plot(x_plot, fit_cubico(x_plot, *popt), '-', color='tab:blue', linewidth=1.5, label='Fit Cubico')
plt.plot(x_plot, f_teorica(x_plot), '--', color='tab:orange', linewidth=1.5, label='Curva Teorica Reticolo')

# Marcatore sull'asse Y per mostrare il limite al continuo ideale (E=1)
plt.axhline(1.0, color='gray', linestyle=':', alpha=0.7)
plt.plot(0, c0, 'ks', label=f'Limite continuo estrapolato ({c0:.4f})')

# Impostazioni assi e legende
plt.xscale('log', base=2)  # Scala logaritmica in base 2 per visualizzare bene i dimezzamenti di a
plt.xlabel(r'$a^2$ (Passo reticolare al quadrato)', fontsize=11)
plt.ylabel(r'$E$ (elemtento matrice efficace)', fontsize=11)
plt.title(r'Estrapolazione al Continuo ($a^2 \to 0$) - Oscillatore Armonico', fontsize=12, pad=15)
plt.legend(loc='lower left', fontsize=10)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)

plt.tight_layout()
plt.show()
