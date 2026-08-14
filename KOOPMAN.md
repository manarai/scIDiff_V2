# Koopman backend — mathematical model and caveats

**scJDO** infers a time-dependent drift field $f_\theta(x,t)$ and its per-cell
Jacobian $J(x,t) = \nabla_x f_\theta$ along pseudotime. The core outputs are
the *temporal Jacobian tensor* $J(\tau) \in \mathbb{R}^{T \times D \times D}$
and its decomposition into interpretable operator archetypes (semi-NMF, default).

The **Koopman backend** is a complementary *spectral lens* on that same
Jacobian sequence. It does **not replace** the drift field, the Jacobian
inference, or scJDO's biological interpretation layer — it adds a linear
operator on the vectorized Jacobian trajectory whose eigenvalues read as
growth / decay rates and oscillation frequencies of the operator dynamics.

The biological semantics still come from scJDO: what changes when Koopman
is enabled is the way the operator sequence is *summarised*, not what
gets measured.

---

## 1. Setup

Let

$$y(\tau) = \operatorname{vec}\bigl(J(\tau)\bigr) \in \mathbb{R}^{D^2}$$

be the flattened Jacobian at pseudotime grid point $\tau$. Stack the
snapshots as rows of a trajectory matrix

$$Y = \bigl[\,y(\tau_1)\;;\;y(\tau_2)\;;\;\dots\;;\;y(\tau_T)\,\bigr]
\in \mathbb{R}^{T \times D^2}.$$

$Y$ is typically wide ($D^2 \gg T$) and low-rank in practice, because the
operator changes smoothly along pseudotime and lives on a small manifold.

---

## 2. Low-rank reduction

Compute the thin SVD

$$Y = U \Sigma V^{\!\top}, \qquad
Z = Y V_r = U_r \Sigma_r \in \mathbb{R}^{T \times r},$$

with $r \ll D^2$. $Z$ is the *reduced* trajectory. All subsequent
Koopman fitting takes place in $\mathbb{R}^r$; the reduction is what
makes EDMD tractable when $D$ is large (a full $D^2 \times D^2$
operator is infeasible for $D \gtrsim 50$).

**Choice of $r$.** Default `reduced_rank = min(T-1, 4·rank+8, D²)`.
Keep $r$ small relative to $T$ so the local least-squares problem is well
posed; enlarging $r$ improves reconstruction but destabilises the local
spectrum.

**Metric caveat.** $Z = Y V_r = U_r \Sigma_r$ carries the singular-value
scaling, so the reduced space has a *non-neutral* metric: directions
with larger singular values dominate the least-squares fit and shape
the eigenvector geometry. `koopman_modes(..., whiten=True)`
orthonormalises the reduced trajectory ($Z = U_r$) so the fit is
balanced across modes. Under exact least squares the two operators are
similar, $K_{\text{whiten}} = \Sigma_r K \Sigma_r^{-1}$, so eigenvalues
coincide and only eigenvectors differ; with `ridge > 0` the regulariser
acts differently in the two metrics so eigenvalues shift slightly too.

Every geometry descriptor in §7 depends on this choice. `eigvec_cond`
obviously does — the eigenvectors themselves change. But `henrici` and
`reactivity` are *also* not invariant under the $\Sigma_r$-similarity:
they are defined against norms ($\lVert \cdot \rVert_F$ and
$\lVert \cdot \rVert_2$), and $\Sigma_r$ is not orthogonal, so a
similarity by $\Sigma_r$ changes the norm even though it does not move
the eigenvalues. `transient_gain` depends on the operator 2-norm and
moves with the metric for the same reason. That is the whole reason
the flag exists. The `whiten` echo in the diagnostics dict ensures the
number and the metric that produced it always travel together.

---

## 3. Extended DMD with identity observables

The Koopman operator $\mathcal{K}$ advances observables, not states:

$$(\mathcal{K}_{\Delta\tau}\, \psi)(y_\tau)
   = \mathbb{E}\bigl[\psi(y_{\tau+\Delta\tau}) \,\big|\, y_\tau\bigr].$$

With identity observables $\psi(y) = y$ this reduces to DMD in the
reduced space. Given snapshot pairs $(z(\tau_i),\, z(\tau_{i+1}))$ the
finite-dimensional approximation solves

$$K \;=\; \arg\min_{K \in \mathbb{R}^{r \times r}}
   \sum_{i} \bigl\| z(\tau_{i+1}) - K\, z(\tau_i) \bigr\|_2^2
   \;+\; \lambda \|K\|_F^2,$$

by ridge-regularised least squares, i.e.
$K = (Z_{-}^{\!\top} Z_{-} + \lambda I)^{-1} Z_{-}^{\!\top} Z_{+}$
where $Z_{-}$ stacks the "before" snapshots and $Z_{+}$ the "after" ones.
The ridge $\lambda$ (parameter `ridge`, default $10^{-4}$) stabilises the
solve when the window is short or nearly rank deficient.

---

## 4. Windowed (nonautonomous) Koopman

Single-cell trajectories are typically **nonstationary** across
pseudotime and branches — a single global $K$ blurs regime-specific
behaviour. We therefore fit a *local* operator per grid point on a
sliding window $W_\tau = [\tau - h,\, \tau + h]$:

$$K_\tau \;=\; \arg\min_{K}
   \sum_{i \in W_\tau} \bigl\| z(\tau_{i+1}) - K\, z(\tau_i) \bigr\|_2^2.$$

A `mode='global'` option fits a single $K$ on all snapshot pairs as a
stationary baseline. The window half-width `window_half` controls the
smoothness / temporal-resolution trade-off; larger windows give a
smoother spectrum at the cost of resolution near sharp transitions.

---

## 5. Eigenvalue interpretation

Eigendecompose each local operator:

$$K_\tau v_k = \lambda_k(\tau)\, v_k.$$

Convert discrete eigenvalues into **continuous-time rates** using the
mean spacing $\Delta\tau$ of that window:

$$\mu_k(\tau) \;=\; \frac{\log \lambda_k(\tau)}{\Delta\tau}
             \;\in\; \mathbb{C}.$$

- $\operatorname{Re}(\mu_k) > 0$ → *growing mode* (local instability of the
  operator sequence at $\tau$).
- $|\operatorname{Im}(\mu_k)| > 0$ → *oscillatory mode* with angular
  frequency $|\operatorname{Im}(\mu_k)|$ radians per unit pseudotime.
- $\operatorname{Re}(\mu_k) < 0$ → contracting / decaying mode.

The dominant modes are those with the largest $|\lambda_k|$. In
`mode='global'` this is the same ordering as the largest
$\operatorname{Re}(\mu_k)$ (single stationary $\Delta\tau$); in
`mode='local'` it is **not** in general — each window has its own
$\Delta\tau$, so $|\lambda|$-ranking and $\operatorname{Re}(\mu)$-ranking
can diverge. Rank on whichever you intend to interpret and be explicit
about the choice. Small-$|\lambda|$ modes tend to alias noise and should
not be over-interpreted (see the branch note below and §8).

**Branch ambiguity of $\log \lambda$.** The complex logarithm is
multi-valued: both $\mu$ and $\mu + 2\pi i k / \Delta\tau$ satisfy
$\lambda = e^{\mu \Delta\tau}$ for every integer $k$. This means the
frequency $|\operatorname{Im}(\mu_k)|$ is *identified* only in
$[0,\, \pi / \Delta\tau)$ radians per unit pseudotime — the **Nyquist
limit**, equivalently $1/(2\Delta\tau)$ cycles per unit pseudotime.
Anything reported above that limit is not measured, it is *aliased*:
principal-branch $\log$ has folded a higher-frequency mode into the
identifiable band. This is the mechanism behind the "large $r$ aliases noise as huge
oscillation frequencies" symptom in §8, but the precise causal chain
is worth spelling out. Small-$|\lambda|$ noise eigenvalues have angles
that are approximately *uniformly distributed* on the unit circle, so
principal-branch $\log \lambda$ gives $|\operatorname{Im}(\log \lambda)|$
values spread across the whole identifiable band $[0, \pi)$ — not
concentrated near the boundary. The observable symptom — implausibly
huge $|\operatorname{Im}(\mu)|$ — comes from the division
$\log \lambda / \Delta\tau$: $|\operatorname{Im}(\log \lambda)| \le \pi$
is bounded, but $\pi / \Delta\tau$ can be arbitrarily large when the
pseudotime step is small, so any noise mode looks like a
high-frequency oscillation once expressed as a continuous-time rate.
The Nyquist companions in the diagnostics dict flag exactly this: if a
reported $|\operatorname{Im}(\mu_k(\tau))|$ sits at the ceiling
`nyquist_ang_freq[τ] = π / Δτ` it is a division artifact, not a
measured oscillation.

`koopman_modes(..., return_diagnostics=True)` returns `delta_tau`,
`nyquist_ang_freq = π / Δτ`, and `nyquist_cycle_freq = 1 / (2Δτ)`
alongside `freqs` so the caveat travels with the numbers. Before
interpreting any oscillation, check that `freqs[t, k]` sits comfortably
below `nyquist_ang_freq[t]`.

Growth $\operatorname{Re}(\mu_k)$ has no such branch ambiguity: only
the imaginary part of $\log \lambda$ is multi-valued. Discretisation
error remains — a strongly growing mode within one $\Delta\tau$ step is
still coarsely sampled — but it is not fundamentally unidentifiable.

---

## 6. Modes and per-time amplitudes

For downstream parity with semi-NMF (patterns × activations), we lift
the reference operator's eigenvectors back to $D \times D$:

$$A_k \;=\; \operatorname{unvec}\bigl(V_r v_k\bigr) \in \mathbb{R}^{D \times D}.$$

Complex conjugate eigenvectors are split into $(\operatorname{Re} v_k,\,
\operatorname{Im} v_k)$ so the returned patterns are real-valued and
consumable by scJDO's existing gene-space projection routines.

Per-time **amplitudes** use the biorthogonal left eigenvectors
$w_k$ (solving $K_\tau^{\!\top} w_k = \lambda_k w_k$ and normalised so
$w_k^{\mathsf{H}} v_k = 1$):

$$b_k(\tau) \;=\; w_k^{\mathsf{H}} z(\tau) \;\in\; \mathbb{C},
   \qquad
   \text{activation}_k(\tau) = |b_k(\tau)|.$$

Magnitudes are non-negative, matching the sign convention downstream
code expects (`instability_scores`, `arch_instability_genes`, …).

**Terminology.** These $A_k$ are Koopman **modes** — eigenvectors of a
linear operator on the reduced trajectory. They are *not* non-negative
regulatory archetypes in the semi-NMF sense, and they carry no built-in
sparsity or biological factoring. The shared `(K, D, D)` return shape
with `jacobian_modes` is an *interface convenience* — it lets the same
downstream code (gene projection, instability scoring) consume either
backend — not a claim that the two objects mean the same thing.
Interpret a semi-NMF pattern as "which program is active"; interpret a
Koopman mode as "which direction the operator is amplifying".

---

## 7. What Koopman gives you

Koopman turns the nonlinear-operator sequence into a **linear spectral
problem on observables**, which makes it especially good at

- extracting slow modes from a mostly-autonomous trajectory,
- separating dominant recurrent structure from transient noise,
- flagging candidate oscillatory components,
- summarising long-horizon behaviour with a small set of rates.

It is *most* useful when there is a clean separation between fast
transients and slower dynamics, and when the pseudotime ordering has
enough resolution to approximate physical snapshot pairs.

### Branch-free geometry descriptors

Alongside the spectral $\lambda / \mu$ picture, the backend returns four
descriptors computed **directly on the discrete-time operator** $K_\tau$
without taking a logarithm. Because they do not touch $\log \lambda$
they are *immune to the branch ambiguity* in §5 — a genuinely
better-posed piece of the output in the snapshot setting than the
frequency estimates. All four depend on the reduced-space metric (see
§2 and the `whiten` flag); document which was used.

| Descriptor | Definition | Reads as |
|---|---|---|
| `henrici`        | $\nu(K) = \sqrt{\lVert K \rVert_F^2 - \sum_i \lvert\lambda_i\rvert^2}\,/\,\lVert K \rVert_F$ | Departure from normality; $\in [0, 1]$, $= 0$ iff $K$ is normal. |
| `reactivity`     | $\sigma_{\max}(K) / \rho(K)$ | One-step amplification vs. asymptotic rate; $\ge 1$, $= 1$ iff normal. |
| `transient_gain` | $\max_{n \in [1, N]} \lVert K^n \rVert_2$ | Short-horizon amplification that asymptotic eigenvalues can miss for non-normal $K$. |
| `eigvec_cond`    | $\kappa(V) = \lVert V \rVert_2 \lVert V^{-1} \rVert_2$ | Warns when the eigen-based interpretation is fragile (nearly defective $K$). |

Large `henrici`, `reactivity`, `transient_gain`, or `eigvec_cond`
values are a signal that the operator is non-normal at that window —
transient amplification, mode mixing, and reactive burst behaviour can
all be present even when every eigenvalue lies inside the unit circle.
Reading only the eigenvalue picture will miss all of it.

**`transient_gain` — two caveats.** Unlike `henrici` and `reactivity`,
which are ratios, `transient_gain` is a norm and inherits two extra
sensitivities the ratios do not have:

- *Horizon-dependent.* It is the maximum of $\lVert K^n \rVert_2$ over
  $n \in [1, N]$ with $N =$ `geometry_n_max` (default 30). Two runs
  with different horizons are not directly comparable — always report
  the horizon alongside the value. If you want a horizon-free summary,
  post-process to the argmax $n^\star$ and store $(n^\star, \text{gain})$
  as a pair.
- *Not scale-free.* Because it is a norm, it moves with any rescaling of
  the reduced trajectory — including the `whiten` metric choice in §2.
  A scale-invariant variant is
  $\lVert K^{n^\star} \rVert_2 \,/\, \rho(K)^{n^\star}$, i.e.
  amplification *relative to the asymptotic rate* at the peak
  horizon; that quantity is invariant under the $\Sigma_r$-similarity
  that whitening induces.

---

## 8. Failure modes and caveats

Koopman is a **spectral summary**, not ground truth. Several caveats
apply specifically in the single-cell setting:

- **Pseudotime is a surrogate.** scRNA-seq is snapshot data — a cell is
  destroyed at measurement. The Koopman operator here summarises the
  *inferred* pseudotime progression, not the intrinsic physical flow.
  Results are more model-dependent and more sensitive to branch
  assignment, time warping, and sampling density than in trajectory
  data.
- **Snapshot data cannot unambiguously identify oscillations.** Treat
  $|\operatorname{Im}(\mu_k)| > 0$ modes as *candidate* oscillations,
  not as physically resolved cycles.
- **Spectral pollution near boundaries.** Windows that clip the
  beginning or end of pseudotime have fewer snapshot pairs and
  systematically noisier spectra. Interior windows are more
  trustworthy; the `peak_interior` helper is a good default guard for
  reported peaks.
- **Not naturally sparse or biologically factored.** Koopman modes are
  mathematically nice but not sparse or non-negative. If your goal is
  interpretable regulatory *programs*, prefer semi-NMF or use both and
  compare.
- **Reduced rank matters.** A too-large $r$ produces spurious
  small-$|\lambda|$ modes that alias noise as huge oscillation
  frequencies; a too-small $r$ misses real spectral structure.

---

## 9. When to use which backend

| Situation | Recommended backend |
|---|---|
| Dense, time-resolved, mostly-autonomous trajectory (differentiation, reprogramming, perturbation recovery) | Koopman + semi-NMF, compared |
| Sparse atlas-like snapshot, strongly branch-heavy fate map | semi-NMF (Koopman spectrum will be noisy) |
| You want to *report a program* (genes, TFs, activation window) | semi-NMF |
| You want to *report a rate* (growth, oscillation, spectral gap) | Koopman |
| Diagnosing whether a transition is one dominant slow mode vs many recurrent switches | Koopman |

---

## 10. API

Three levels — pick the one that matches how deep you want to go.

**High-level (main API).** Same `fit_drift` call, switch the backend:

```python
sjd.tl.fit_drift(
    adata,
    archetype_method="koopman",
    koopman_kwargs=dict(mode="local", window_half=8, ridge=1e-4),
)
res  = adata.uns["scjdo"]
koop = res["koopman"]                     # spectral diagnostics
koop["eigenvalues"]                       # (T, K) complex λ_k(τ)
koop["growth_rates"]                      # (T, K) Re(log λ / Δτ)
koop["freqs"]                             # (T, K) |Im(log λ / Δτ)| — check against nyquist
koop["delta_tau"], koop["nyquist_ang_freq"], koop["nyquist_cycle_freq"]
koop["geometry"]                          # branch-free: henrici, reactivity,
                                          # transient_gain, eigvec_cond, each (T,)
koop["whiten"]                            # metric echo (see §2)
```

`res["patterns"]`, `res["activations"]`, `res["max_real_eig"]` and every
existing downstream field (gene scores, instability tables, regulator
inference) keep working — the Koopman backend returns the same
`(K × D × D, T × K, error)` shape as semi-NMF by construction.

**Dispatcher on a pre-computed tensor.** When you already hold a
Jacobian tensor and just want to compare backends:

```python
patterns, activations, err, diag = sjd.tl.decompose_archetypes(
    J_tensor, method="koopman", rank=5,
)
```

**Low-level.** Direct access to all Koopman options:

```python
from scjdo.archetypes import koopman_modes

patterns, activations, err, diag = koopman_modes(
    J_tensor,
    rank=5,
    reduced_rank=None,       # auto: min(T-1, 4·rank+8, D²)
    mode="local",            # or "global"
    window_half=8,
    ridge=1e-4,
    grid=t_centers,          # enables continuous-time rate conversion
    smooth_sigma=1.0,
    whiten=False,            # metric choice; see §2
    geometry_n_max=30,       # transient-gain horizon
    return_diagnostics=True,
)
```

---

## 11. Reproducibility

The Koopman backend is deterministic given the reduced trajectory
(SVD + least squares are both deterministic to fp precision). Randomness
enters only via the *upstream* drift training and the kernel-windowing
bootstrap; both are seeded and preserved. Store `adata.uns[key]["koopman"]`
alongside the semi-NMF outputs for reproducible comparison.

---

## 12. See also

- [`MATH.md`](MATH.md) — full derivation of the drift field, Jacobian
  tensor, kernel windowing, and semi-NMF archetype decomposition.
- `scjdo.archetypes.koopman` — implementation, docstring, and inline
  references to the equations above.
- `tests/test_koopman.py` — shape, spectral-recovery, and
  end-to-end integration tests.
