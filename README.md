# Quantum Harmonic Oscillator on the Lattice

A C implementation of the quantum harmonic oscillator using lattice field theory techniques: path integral formulation, Metropolis Monte Carlo, and jackknife error estimation.

---

## Table of Contents

- [Physics Background](#physics-background)
- [Numerical Methods](#numerical-methods)
- [Project Structure](#project-structure)
- [Modules](#modules)
- [Parameters](#parameters)
- [How to Build and Run](#how-to-build-and-run)
- [Output Files](#output-files)
- [Continuum Limit](#continuum-limit)
- [Dependencies](#dependencies)

---

## Physics Background

The quantum harmonic oscillator is one of the few systems with an exact analytical solution, making it an ideal benchmark for lattice field theory methods.

### Euclidean Action

Working in imaginary time (Euclidean formulation) with ℏ = 1, the continuum action is:

```
S_E[x] = ∫_0^T dτ [ (M/2) ẋ² + (M ω²/2) x² ]
```

with periodic boundary conditions x(0) = x(T).

### Lattice Discretisation

The time interval [0, T] is divided into N steps of size a, so that T = N·a. The path integral becomes a finite-dimensional integral over the lattice variables x_i ≡ x(i·a):

```
S_lat = Σ_{i=0}^{N-1} [ (M/2a)(x_{i+1} - x_i)² + (a M ω²/2) x_i² ]
```

with the periodic condition x_N ≡ x_0. In terms of the adimensional (lattice) variables M̂ = m/a and Ω = ωa, this becomes:

```
S_lat = Σ_{i=0}^{N-1} [ (M̂/2)(x_{i+1} - x_i)² + (M̂ Ω²/2) x_i² ]
```

### Observables

The two-point correlation function of the position operator is:

```
G(t) = (1/N) Σ_{k=0}^{N-1} x_k · x_{(k+t) mod N}
```

From G(t) one extracts the energy gap and the matrix element:

**Energy gap** (effective mass):

```
a·ΔE(t) = arccosh[ (G(t+1) + G(t-1)) / (2 G(t)) ]
```

**Matrix element** of the position operator between ground and first excited state:

```
|⟨E_0|x̂|E_1⟩|²(t) = G(t) · exp(N/2 · a·ΔE(t)) / [ 2 cosh((N/2 - t) · a·ΔE(t)) ]
```

### Continuum Limit

The Taylor expansion of the observables in powers of a² gives:

```
ΔE(a) = ω [ 1 - (ωa)²/24 + 3(ωa)⁴/640 - 5(ωa)⁶/7168 + O(a⁸) ]

|⟨E_0|x̂|E_1⟩|(a) = 1/√(2mω) [ 1 - (ωa)²/16 + 5(ωa)⁴/512 + O(a⁶) ]
```

The continuum value is recovered by fitting a polynomial in a² and extrapolating to a = 0.

---

## Numerical Methods

### Metropolis Algorithm

The Metropolis algorithm generates configurations distributed according to the Boltzmann weight e^{-S}. Each sweep updates all N lattice sites in sequence:

```
for each site i = 0, ..., N-1:
    propose  x'_i = x_i + 2Δ·(r₁ - 0.5)      r₁ ~ Uniform[0,1)
    compute  ΔS = S(x') - S(x)
    draw     r₂ ~ Uniform[0,1)
    if exp(-ΔS) > r₂:
        accept  x_i ← x'_i
```

The step size Δ is chosen so that the acceptance rate is approximately 50%. Following Creutz & Freedman, the optimal value is Δ = 2√a.

The variation of the action when only site i changes (used for efficiency) is:

```
ΔS = (M̂/a) δx [ (2 + Ω²)(x_i + δx/2) - x_{i-1} - x_{i+1} ]
```

where δx = x'_i - x_i.

### Binning

Successive Monte Carlo configurations are correlated. The autocorrelation time τ is estimated by fitting the normalised autocorrelation function:

```
Γ(t_M) = [ ⟨C_n(t) C_{n+t_M}(t)⟩ - ⟨C_n(t)⟩² ] / [ ⟨C_n²(t)⟩ - ⟨C_n(t)⟩² ]
```

The bin size is chosen as Bin ≈ 10·τ so that binned measurements are statistically independent. With n = N_conf / Bin independent bins, the effective number of independent measurements is:

```
n_eff = N_conf / (2τ)
```

### Jackknife Error Estimation

Because ΔE is a non-linear function of G(t), errors cannot be propagated linearly. The jackknife method is used instead.

Given n binned measurements G^(0)(t), ..., G^(n-1)(t):

1. **Global mean**: G̅(t) = (1/n) Σ_k G^(k)(t)

2. **Jackknife mean** (leaving out bin k):
```
G̅^(k)(t) = [ n·G̅(t) - G^(k)(t) ] / (n - 1)
```

3. **Observable on global mean**:
```
ΔE̅(t) = arccosh[ (G̅(t+1) + G̅(t-1)) / (2 G̅(t)) ] / a
```

4. **Observable on jackknife sample**:
```
ΔE^(k)(t) = arccosh[ (G̅^(k)(t+1) + G̅^(k)(t-1)) / (2 G̅^(k)(t)) ] / a
```

5. **Jackknife error**:
```
σ_{ΔE}(t) = √[ (n-1)/n · Σ_k (ΔE^(k)(t) - ΔE̅(t))² ]
```

### Weighted Average

The final value of each observable is obtained as a weighted average over the plateau region (points with relative error ≤ 10%):

```
ΔE_final = [ Σ_t ΔE(t)/σ²(t) ] / [ Σ_t 1/σ²(t) ]

σ_final = 1 / √[ Σ_t 1/σ²(t) ]
```

---

## Project Structure

```
LFC/
├── include/              # Header files
│   ├── global.h          # Physical parameters and global arrays
│   ├── action.h          # Function declarations
│   ├── random.h          # RANLUX random number generator
│   ├── start.h           # Initialisation routines
│   └── types.h           # Type definitions
│
├── modules/
│   ├── action/           # Physics modules
│   │   ├── soa.c                      # Total action S[x]
│   │   ├── ds.c                       # Action variation ΔS
│   │   ├── metropolis.c               # Metropolis sweep
│   │   ├── correlatore.c              # Two-point correlator G(t)
│   │   ├── autocorrelatore.c          # Autocorrelation function Γ(t_M)
│   │   ├── correlatore_mean.c         # Mean of binned correlator
│   │   ├── correlatore_sigma.c        # Std dev of binned correlator
│   │   ├── effective_mass_jackknife.c # ΔE(t) with jackknife errors
│   │   ├── energy.c                   # ΔE from global mean
│   │   ├── energy_jackknife.c         # Jackknife error on ΔE
│   │   └── matrix_element.c           # |⟨E_0|x̂|E_1⟩|²(t)
│   ├── random/           # RANLUX random number generators
│   │   ├── ranlxs.c      # Single precision
│   │   ├── ranlxd.c      # Double precision
│   │   └── gauss.c       # Gaussian random numbers
│   └── start/            # Utility routines
│       ├── start.c        # Generator initialisation
│       └── utils.c        # Memory allocation, error handling
│
├── main/
│   ├── main.c            # Main simulation program
│   └── Makefile
│
├── devel/
│   ├── oa/               # Tests for action modules
│   │   ├── test_action.c     # Validates soa.c
│   │   ├── test_ds.c         # Validates ds.c
│   │   ├── test_acceptance.c # Measures Metropolis acceptance rate
│   │   └── test_plateau.c    # Analyses the effective mass plateau
│   └── random/           # Tests for random number generators
│
├── plots/                # Python analysis scripts
│   ├── plot_action.py         # Thermalisation plot
│   ├── plot_autocorr.py       # Autocorrelation comparison
│   └── limite_continuo.py     # Continuum limit fit
│
└── doc/                  # Documentation
```

---

## Modules

### `soa.c` — Total Action

Computes the full lattice action:

```
S = Σ_{i=0}^{N-1} [ (M̂/2)(x_{i+1 mod N} - x_i)² + (M̂ Ω²/2) x_i² ]
```

```c
double soa(void);
```

### `ds.c` — Local Action Variation

Computes the change in action when site i is updated by δx, without recomputing the full action. Only three terms depend on x_i:

```
ΔS = M̂ δx [ (2 + Ω²)(x_i + δx/2) - x_{i-1} - x_{i+1} ]
```

```c
double ds(int i, double dxx);
```

### `metropolis.c` — Metropolis Sweep

Performs one full sweep over all N lattice sites using the Metropolis algorithm with step size DELTA.

```c
void metropolis(void);
```

### `correlatore.c` — Two-Point Correlator

Computes the normalised two-point function for all time separations t:

```
G(t) = (1/N) Σ_{k=0}^{N-1} x_k · x_{(k+t) mod N}
```

Returns a dynamically allocated array of N doubles (caller must free).

```c
double *correlatore(void);
```

### `autocorrelatore.c` — Autocorrelation Function

Computes the normalised autocorrelation Γ(t_M) of a time series up to lag tmax and writes the result to a text file.

```c
void autocorrelatore(double corr[], int n, int tmax, char *filename);
```

### `correlatore_mean.c` — Binned Correlator Mean

Reads the binary file of binned correlators and returns the mean G̅(t) over all bins.

```c
double *correlatore_mean(const char *filename_bin);
```

### `correlatore_sigma.c` — Binned Correlator Standard Deviation

Reads the binary file of binned correlators and returns the standard deviation of the mean σ(t).

```c
double *correlatore_sigma(const char *filename_bin, double *media);
```

### `effective_mass_jackknife.c` — Effective Mass with Jackknife Errors

Computes ΔE(t) from the global mean and its jackknife error σ_{ΔE}(t) by reading the binned correlator file one row at a time.

```c
void effective_mass_jackknife(double *corr_mean, double *dE, double *sigma_dE);
```

### `matrix_element.c` — Matrix Element

Computes |⟨E_0|x̂|E_1⟩|²(t) from the mean correlator and the effective mass:

```
M(t) = G̅(t) · exp(ΔE(t)·N·a/2) / [ 2 cosh(ΔE(t)·a·(N/2 - t)) ]
```

```c
double *matrix_element(double *corr_mean, double *dE);
```

---

## Parameters

All physical and simulation parameters are defined in `include/global.h`:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `N` | N | Number of lattice sites |
| `a` | a | Lattice spacing |
| `M` | m | Particle mass (physical) |
| `W` | ω | Oscillator frequency (physical) |
| `DELTA` | Δ | Metropolis step size (≈ 2√a) |
| `Nconf` | N_conf | Number of Monte Carlo sweeps |
| `Bin` | Bin | Bin size (≈ 10·τ) |
| `Nacor` | N_acor | Maximum lag for autocorrelation |

The lattice parameters are:

```
M̂ = m/a      (lattice mass)
Ω  = ω·a     (lattice frequency)
T  = N·a     (total Euclidean time, kept fixed in continuum limit)
```

### Simulation parameters used for the continuum limit (T = 64)

| N | a | Δ | Bin | τ | N_conf |
|---|---|---|-----|---|--------|
| 32 | 2.000 | 1.60 | 25 | 1.3 | 125 000 |
| 64 | 1.000 | 1.70 | 50 | 2.6 | 250 000 |
| 128 | 0.500 | 1.20 | 100 | 8.2 | 500 000 |
| 256 | 0.250 | 1.00 | 300 | 31 | 1 500 000 |
| 512 | 0.125 | 0.75 | 1000 | 93 | 5 000 000 |
| 1024 | 0.0625 | 0.50 | 5000 | 514 | 25 000 000 |

---

## How to Build and Run

### Requirements

- GCC with C99 support
- GNU Make
- Python 3 with NumPy, SciPy, Matplotlib (for analysis scripts)

### Build

```bash
cd main/
make
```

To clean:

```bash
make clean
```

### Run

```bash
./main
```

The program:
1. Thermalises the system (500 sweeps from a cold start x_i = 0)
2. Runs N_conf sweeps, computing and binning the correlator
3. Writes the binned correlator to `correlatorebinn.bin`
4. Computes autocorrelation functions
5. Computes the effective mass ΔE(t) and its jackknife error
6. Computes the matrix element |⟨E_0|x̂|E_1⟩|²(t)
7. Outputs a weighted average for each observable

### Running the tests

```bash
cd devel/oa/
make
./test_action       # validates soa.c against analytical values
./test_ds           # validates ds.c against finite differences
./test_acceptance   # measures the Metropolis acceptance rate
./test_plateau      # analyses the effective mass plateau
```

### Analysis

```bash
cd plots/
python3 plot_action.py        # thermalisation plot
python3 plot_autocorr.py      # autocorrelation comparison (binned vs raw)
python3 limite_continuo.py    # continuum limit fit
```

---

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `termalizzazione.txt` | ASCII: sweep S | Action vs sweep number |
| `correlatorebinn.bin` | Binary: float64 | Binned correlator matrix [n_bin × N] |
| `medie_sigma_corr.txt` | ASCII: t G̅(t) σ(t) | Mean and std dev of correlator |
| `autocorrelatore.txt` | ASCII: t_M Γ(t_M) | Autocorrelation of raw series |
| `autocorrelatorebin.txt` | ASCII: t_M Γ(t_M) | Autocorrelation of binned series |
| `effective_mass.txt` | ASCII: t ΔE(t) σ(t) | Effective mass with jackknife errors |

---

## Continuum Limit

The continuum limit is extracted by simulating at several values of a (keeping T = N·a fixed) and fitting the observables as polynomials in a²:

```
ΔE(a²) = c_0 + c_1·a² + c_2·a⁴ + c_3·a⁶
```

The coefficient c_0 gives the continuum value. The analytical predictions are:

```
c_0 = ω        = 1.0000
c_1 = -ω³/24  = -0.04167
c_2 = 3ω⁵/640 = +0.00469
c_3 = -5ω⁷/7168 = -0.00070
```

The fit is performed with `plots/limite_continuo.py` using `scipy.optimize.curve_fit` with `absolute_sigma=True`.

---

## Dependencies

- **RANLUX** random number generator (M. Lüscher): included in `modules/random/`
- **Python packages**: `numpy`, `scipy`, `matplotlib`

---

## References

1. M. Creutz and B. Freedman, *A statistical approach to quantum mechanics*, Ann. Phys. 132 (1981) 427.
2. H. J. Rothe, *Lattice Gauge Theories: An Introduction*, World Scientific (2012).
3. M. Lüscher, *A portable high-quality random number generator for lattice field theory simulations*, Comput. Phys. Commun. 79 (1994) 100.
