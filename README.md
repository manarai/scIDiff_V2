# 🧬 scIDiff: Schrödinger Bridge Learning of Single-Cell Regulatory Dynamics

scIDiff (single-cell inverse Diffusion) learns time-dependent gene-regulatory drift fields directly from single-cell data. It combines Schrödinger Bridge optimal transport, score-based generative modeling, and regulatory network inference to explain how cells move through gene-expression space over time.

## 🌌 Mathematical Foundations

Given single-cell measurements across developmental or perturbation conditions, scQDiff models a continuous stochastic process:

$$
dX_t = u(X_t, t)dt + \sqrt{2\beta}dW_t
$$

- **$X_t \in \mathbb{R}^d$**: cell state (e.g., gene expression)
- **$u(X_t,t)$**: drift field (deterministic regulatory flow learned by scQDiff)
- **$\beta$**: diffusion constant  
- **$dW_t$**: Brownian noise

The drift is learned so simulated trajectories transport the empirical initial distribution $\rho_0$ to the terminal distribution $\rho_1$ while minimizing a Schrödinger Bridge objective:

$$
\min_u \; \text{KL}(P_u \| P_{\text{ref}}) \quad \text{subject to} \quad \rho(0) = \rho_0, \; \rho(1) = \rho_1.
$$

A useful summary is the temporal Jacobian:

$$
J(t) = \frac{\partial u}{\partial x}(t)
$$

which encodes local (causal) gene-to-gene influence along time.

The model also exposes the expected state drift:

$$
u(x,t) \approx \frac{d\mathbb{E}[X_t]}{dt}.
$$

## 🔍 What scQDiff Learns

| Layer | Description | Output |
|-------|-------------|--------|
| Regulatory drift | Deterministic direction of transcriptional change | $u(x,t)$ |
| Temporal Jacobian | Local gene→gene influence (causal) | $J(t)=\tfrac{\partial u}{\partial x}$ |
| Regulatory archetypes | Low-rank temporal modes of regulation | $u(x,t)\approx\sum_k a_k(t)A_k x$ |
| (Optional) Communication | Dynamic ligand–receptor signaling (CellPhoneDB) | $W_{ij}(t)$ |

## 🔗 Extension: scQDiff-Comm (Cell–Cell Communication)

With CellPhoneDB ligand–receptor priors, scQDiff can model time-evolving communication graphs:

$$
dX_t^{(i)} = u_{\text{intra}}(X_t^{(i)}, t)dt + \sum_j W_{ij}(t) f(X_t^{(j)} - X_t^{(i)})dt + \sqrt{2\beta}dW_t^{(i)}.
$$

- **$W_{ij}(t)$**: communication strength between cells $i$ and $j$ at time $t$
- **$f(\Delta x)$**: interaction function (e.g., linear or gated linear)

Result: communication archetypes (e.g., inflammatory relay, exhaustion/resolution) over time.

## 🧮 Core Features

- Time-continuous modeling of differentiation/perturbation trajectories
- Regulatory Jacobians for directional, causal interpretation  
- Forward & reverse simulation for counterfactual/reprogramming paths
- Cell–cell communication dynamics via CellPhoneDB (optional)
- Cytoscape exports for dynamic networks (.graphml, .cyjs)
- Multi-omics ready (ATAC, RNA velocity, protein/metabolite embeddings)

## 💻 Installation

```bash
# Create conda env
conda create -n scqdiff python=3.10
conda activate scqdiff

# Clone the repository
git clone https://github.com/manarai/scQDiff.git
cd scQDiff

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
