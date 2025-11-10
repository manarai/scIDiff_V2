# 🧬 scIDiff: Schrödinger Bridge Learning of Single-Cell Regulatory Dynamics

**scIDiff** (*single-cell inverse Diffusion*) learns **time-dependent gene-regulatory drift fields** from single-cell data by solving a **Schrödinger Bridge (SB)** problem.  
It unifies **optimal transport**, **score-based generative modeling**, and **regulatory network inference** to reconstruct how cells **flow through gene-expression space** over time — and how **irreversibility** emerges from these trajectories.

---

## 🌌 Mathematical Foundations

Given observed cell states sampled at different times or perturbation conditions, scIDiff models a **stochastic process**:

\$$
dX_t = u(X_t, t)\,dt + \sqrt{2\beta}\,dW_t
\$$
- **$X_t \in \mathbb{R}^d$** — cell state (e.g., transcriptome)
- **$u(X_t,t)$** — drift field (regulatory flow learned by scIDiff)
- **$\beta$** — diffusion constant  
- **$dW_t$** — Brownian noise  

The learned drift transports the empirical initial distribution $\rho_0$ to the terminal $\rho_1$ while minimizing a Schrödinger-Bridge energy:

\$$
\min_u \; \text{KL}(P_u \,\|\, P_{\text{ref}}) \quad \text{s.t.} \quad \rho(0)=\rho_0, \;\rho(1)=\rho_1.
\$$

A **temporal Jacobian**

\$$
J(t) = \frac{\partial u}{\partial x}(t)
\$$

encodes local, causal gene-to-gene influence along time, revealing how gene modules control differentiation flow.

---

## 🔁 Forward and Reverse Drift

scIDiff learns both **forward** and **reverse** drift fields:

$$
\begin{aligned}
\text{Forward:} \quad & dX_t = u(x,t)\,dt + \sqrt{2\beta}\,dW_t \\
\text{Reverse:} \quad & dX_t = \big[u(x,t) - 2\beta\nabla_x \log\rho_t(x)\big]\,dt + \sqrt{2\beta}\,d\bar{W}_t
\end{aligned}
$$

- The **forward drift** $u(x,t)$ describes how cells evolve naturally.  
- The **reverse drift** captures how much “work” would be required to **reprogram** cells backward in time (e.g., iPSC or rejuvenation).  
- The **difference** $\Delta u = u_{\text{fwd}} - u_{\text{rev}}$ quantifies **irreversibility** — an analog of **biological entropy production**.

This forward–reverse asymmetry provides a rigorous way to identify **irreversible cell-fate decisions** and the regulators (Yamanaka-like factors) capable of reversing them.

---

## 🧠 What scIDiff Learns

| Layer | Description | Output |
|:------|:-------------|:--------|
| **Regulatory drift** | Deterministic direction of transcriptional change | $u(x,t)$ |
| **Reverse drift** | Reprogramming direction opposing natural evolution | $u_{\text{rev}}(x,t)$ |
| **Irreversibility field** | Entropic asymmetry between forward and reverse | $\Delta u = u - u_{\text{rev}}$ |
| **Temporal Jacobian** | Local causal gene-to-gene influence | $J(t)=\tfrac{\partial u}{\partial x}$ |
| **Regulatory archetypes** | Low-rank modes of temporal regulation | $u(x,t)\approx\sum_k a_k(t)A_kx$ |
| **(Optional) Communication** | Dynamic ligand–receptor signaling between cells | $W_{ij}(t)$ |

---

## 🔗 Extension: scIDiff-Comm (Cell–Cell Communication)

With **CellPhoneDB** ligand–receptor priors, scIDiff models **time-evolving communication graphs**:

$$
dX_t^{(i)} = u_{\text{intra}}(X_t^{(i)}, t)dt
+ \sum_j W_{ij}(t)f(X_t^{(j)} - X_t^{(i)})dt
+ \sqrt{2\beta}dW_t^{(i)}.
$$

- **$W_{ij}(t)$** — communication strength between cells  
- **$f(\Delta x)$** — interaction kernel  

Result: **communication archetypes** (e.g., inflammatory relay, exhaustion/resolution) that co-evolve with regulatory drift.

---

## ⚖️ Quantifying Irreversibility

scIDiff computes **entropy production** and **cycle flux** metrics to measure the degree of biological irreversibility:

\$$
\dot{S}(t) = \mathbb{E}\!\left[\frac{\|u(x,t)-u_{\text{rev}}(x,t)\|^2}{2\beta}\right]
\$$

- High $\dot{S}(t)$ → irreversible differentiation (e.g., commitment, exhaustion)  
- Low $\dot{S}(t)$ → reversible or plastic states (e.g., stem, progenitor)  

This enables principled identification of **control nodes** capable of restoring reversibility.

---

## 🧮 Core Features

- Schrödinger-Bridge learning of forward & reverse dynamics  
- Directional Jacobians for causal inference  
- Quantitative irreversibility (entropy production, Δu)  
- Counterfactual simulations for reprogramming paths  
- Multi-omics ready (RNA, ATAC, velocity, protein, metabolite embeddings)  
- Cytoscape exports for dynamic networks (`.graphml`, `.cyjs`)  
- Optional scIDiff-Comm module for cell–cell signaling

---

## 💻 Installation

```bash
# Create conda environment
conda create -n scidiff python=3.10
conda activate scidiff

# Clone repository
git clone https://github.com/manarai/scIDiff.git
cd scIDiff

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
