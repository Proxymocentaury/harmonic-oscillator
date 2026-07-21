# Lattice Feynman Correlator

Simulazione Monte Carlo (path integral Euclideo) dell'oscillatore armonico quantistico su reticolo, per l'estrazione dello spettro energetico (gap fondamentale) e dell'elemento di matrice della posizione, con successiva estrapolazione al limite continuo.

## Fisica del problema

Si discretizza il path integral Euclideo di un oscillatore armonico 1D su un reticolo temporale periodico di `N` siti e passo reticolare `a`, con azione

```
S = Σ_i [ M/(2a)·(x_i − x_{i−1})² + (a·M·W²/2)·x_i² ]
```

Campionando le configurazioni `{x_i}` con un algoritmo Metropolis locale si stima il correlatore a due punti

```
G(t) = ⟨x(0) x(t)⟩ = (1/N) Σ_k x_k · x_{(k+t) mod N}
```

da cui si estraggono, tramite jackknife:

- **ΔE = E₁ − E₀**, il gap tra stato fondamentale e primo eccitato, tramite la massa efficace locale `acosh[(G(t+1)+G(t−1))/(2G(t))] / a` (esatta per una forma a singolo esponenziale periodizzata, non solo approssimata);
- **|⟨0|x|1⟩|²**, l'elemento di matrice, dal fit della forma funzionale periodica `G(t) = M·[e^{−E t} + e^{−E(N−t)}]`.

Il limite al continuo (`a → 0` a `T = N·a` fissato) va confrontato con le predizioni esatte non perturbative

```
cosh(E·a) = 1 + a²W²/2        Mth = a / (2·sinh(E·a))
```
che si riducono a `E → W` e `M → 1/(2MW)` per `a → 0`, e con lo sviluppo perturbativo in potenze di `a²` usato negli script di fit (valido solo per `a` sufficientemente piccolo).

## Struttura del progetto

```
LFC/
├── include/            # header condivisi (parametri globali, prototipi)
│   ├── global.h         # N, a, M, W, DELTA, Nconf, Bin, Nacor (compile-time)
│   ├── action.h
│   ├── random.h
│   └── start.h
├── modules/
│   ├── action/          # azione, Metropolis, correlatori, jackknife
│   ├── random/          # generatore RANLUX (ranlxs/ranlxd) + gaussiano
│   └── start/           # (non utilizzato dal main attuale)
├── main/
│   ├── main.c           # loop di termalizzazione + Monte Carlo + analisi
│   └── Makefile
├── plots/               # script Python per analisi/estrapolazione
│   ├── plot_action.py          # verifica termalizzazione
│   ├── plot_autocorr.py        # stima del tempo di autocorrelazione τ
│   ├── correloratore_stat.py
│   ├── limite_continuo.py      # fit ΔE(a²) → continuo
│   └── M_limite_continuo.py    # fit elemento di matrice(a²) → continuo
└── doc/                 # documentazione RANLUX
```

## Compilazione

```bash
cd main
make
```

Il Makefile genera automaticamente le dipendenze (`.d`), quindi modificare `include/global.h` forza la ricompilazione di tutti i moduli che ne dipendono al `make` successivo.

## Parametri della simulazione

Tutti i parametri fisici e statistici sono **costanti di compilazione** in `include/global.h`:

| Parametro | Significato |
|---|---|
| `N` | numero di siti reticolari (estensione Euclidea in unità di `a`) |
| `a` | passo reticolare |
| `M`, `W` | massa e frequenza dell'oscillatore (unità naturali) |
| `DELTA` | ampiezza della proposta Metropolis (tarata per accettazione ≈ 50%) |
| `Nconf` | numero totale di sweep di misura |
| `Bin` | dimensione del bin per il binning/jackknife (≈ 10·τ) |
| `Nacor` | numero di lag per cui si calcola la funzione di autocorrelazione |

**Nota importante**: `N·a` deve restare **fissato** su tutti i lanci (rappresenta l'estensione fisica Euclidea T), variando solo `N` e `a` in modo inversamente proporzionale, per poter fare un confronto valido nel limite al continuo.

## Workflow per il limite al continuo

1. **Stima τ**: per ogni valore di `a` che si vuole includere nella scansione, fai un run preliminare con `Nacor` grande (es. 100–500) e `Nconf` moderato, poi stima il tempo di autocorrelazione integrato con `plot_autocorr.py` (fit esponenziale della funzione di autocorrelazione).
2. **Imposta `Bin ≈ 10·τ`** (di conseguenza `NTHERM = 20·Bin`) e scegli `Nconf` in modo che `Nconf/Bin` (il numero di misure effettivamente indipendenti) sia lo stesso per ogni reticolo.
3. **Tara `DELTA`** per un tasso di accettazione vicino al 50%.
4. Ricompila (`make`) e lancia (`./main`) per ciascuna coppia `(N, a)` con `N·a` costante.
5. Verifica `plots/termalizzazione.png` e `plots/autocorrelatore.png` (decadimento pulito).
6. Colleziona `ΔE(a)` e `M(a)` da tutti i run in `plots/dati.txt`, poi lancia `limite_continuo.py` / `M_limite_continuo.py` per il fit polinomiale in `a²` e l'estrapolazione a `a → 0`.

## Output prodotti da un run

| File | Contenuto |
|---|---|
| `termalizzazione.txt` / `.png` | azione durante la fase di termalizzazione |
| `autocorrelatore.txt`, `autocorrelatorebin.txt` | funzione di autocorrelazione normalizzata, sulle misure grezze e sui bin |
| `correlatorebinn.bin` | correlatore binnato, formato binario (`double[N]` per bin) — usato dal jackknife |
| `medie_sigma_corr.txt` | media e deviazione standard naive del correlatore (controllo incrociato) |
| `effective_mass.txt` | massa efficace `dE[t]` ± errore jackknife, per ogni `t` |
| `matrix_element.txt` | elemento di matrice `mat_el[t]` ± errore jackknife |
| stdout | `ΔE` e `Mel` finali (media pesata sul plateau con taglio su errore relativo ≤ 10%) |

## Limitazioni note / lavori in corso

- Il calcolo di `correlatore()` è O(N²) per sweep; una versione basata su FFT (tutti gli `N` usati sono potenze di 2) ridurrebbe il costo a O(N log N).
- La combinazione finale dei punti di plateau `dE[t]`/`mat_el[t]` in `main.c` usa una media pesata assumendo indipendenza tra i vari `t`, che in realtà sono correlati (stessi bin jackknife) — l'errore finale su `ΔE`/`Mel` è probabilmente leggermente sottostimato rispetto a un jackknife applicato all'intera combinazione.
- normalizzazione del passo reticolare

## Dipendenze

- Compilatore C compatibile con C99 (`gcc`)
- Python 3 con `numpy`, `scipy`, `matplotlib` per gli script di analisi in `plots/`
