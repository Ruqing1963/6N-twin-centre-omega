# 6N Twin-Centre Factor Statistics (Part XII)

What distribution does omega₍>3₎(N) follow for the sparse set of twin centres
themselves? A single clean replacement: **1/q → 1/(q-2)**.

**The single-prime probability.** An ordinary integer is divisible by a prime q
with probability 1/q. A twin centre N must avoid the two forbidden residues
dead(q) = {±6⁻¹ mod q} (so 6N∓1 can both be prime), leaving q-2 admissible
residues — of which 0 (i.e. q|N) is one. Hence:

```
    P(q | N , N twin) = 1/(q-2)        (q > 3)
```

versus 1/q for ordinary integers. Verified per prime on S₉ and S₁₀
(23,988,173 twin centres): for small q the measured probability equals 1/(q-2) to
within 0.4% (0.1% for the smallest, most frequent q), tightening from S₉ to S₁₀.
Large q are sampling-limited (P~1/q gives only tens of divisible centres at
q~10⁵), not a deviation from the law.

**Root of the Part I enrichment.** The boosted probability gives the relative
density factor

```
    [P(q|N,twin)/P(q∤N,twin)] / [(1/q)/(1-1/q)] = (q-1)/(q-3),
```

so the Part I enrichment ∏(q-1)/(q-3) is exactly the product of these single-prime
probability boosts 1/q → 1/(q-2). The macroscopic enrichment is the microscopic
divisibility law, multiplied up.

**Mean shift.** Since omega₍>3₎(N) = Σ_{q>3} 1{q|N},

```
    <omega>_twin - <omega>_ord = Σ_{q>3} (1/(q-2) - 1/q) = 0.2604,
```

a convergent constant (each term ~2/q²), cutoff-independent — unlike the
individual means, which diverge as lnln. The measured finite-N shift grows with
the shell (+0.2220 on S₉, +0.2270 on S₁₀), approaching the closed form as the
large-q truncation recedes.

> **Normality is NOT claimed.** At the accessible scale lnln(6N)≈3 the
> distribution has not reached its Erdős–Kac limit: mean (2.70) ≠ variance (1.04),
> skew +0.26, excess kurtosis −0.29 — for ordinary integers as much as for twins.
> Whether omega₍>3₎ over twin centres tends to a normal law with Bernoulli
> parameter 1/(q-2), and what governs its (constrained) variance, is left as an
> **open problem**.
>
> **Scope.** Experimental / computational number theory; S₉, S₁₀. No claim about
> the infinitude of twin primes or any k-tuple conjecture.

Part I: doi:10.5281/zenodo.20470367 · VIII: doi:10.5281/zenodo.20519998 ·
XI: doi:10.5281/zenodo.20526835

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   ├── omega_perq_data.csv  q, P_ord, P_tw_S9, P_tw_S10, inv_qm2, ratio_S10
│   ├── omega_dist_data.csv  k, ord_S9, twin_S9  (omega distribution histogram)
│   └── omega_summary.csv    per-shell mean/var/skew/kurtosis and mean shift
├── code/
│   ├── omega_twin.py        per-prime P(q|N|twin) vs 1/(q-2), mean shift vs closed
│   │                        form, normality diagnostics; emits omega_perq/summary CSV
│   └── make_omega_fig.py    builds the 2-panel figure from ../data
├── figures/                fig_paper12_omega.{pdf,png}
└── paper/                  Chen_6N_Paper12.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Per-prime probabilities, mean shift, normality diagnostics. Default S10
#    (~17 min). Emits omega_perq_S{K}.csv and omega_summary_S{K}.csv.
python code/omega_twin.py            # S10
MAXK=9 python code/omega_twin.py     # S9 (shows the mean shift is smaller)

# 2. Figure (reads ../data/omega_perq_data.csv and omega_dist_data.csv).
cd code && python make_omega_fig.py
```

### Conventions (same as Parts I–XI)

- Twin centre N: 6N−1, 6N+1 both prime. omega₍>3₎(N) = # distinct prime factors >3.
- dead(q) = {±6⁻¹ mod q}: residues making 6N∓1 divisible by q. q-2 admissible residues.
- P(q|N|twin): fraction of twin centres divisible by q (measured); equals 1/(q-2).
- Mean shift closed form Σ_{q>3}(1/(q-2)−1/q) = 0.2604 (convergent).
- Engine: complete segmented-sieve factorisation (unbiased omega) + interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
