# 🧬 scIDiff: Schrödinger Bridge Learning of Single-Cell Regulatory Dynamics

How do cells decide their fate? What are the key genes that drive differentiation, and can we reverse the process? `scIDiff` (single-cell Inverse Diffusion) is a powerful framework designed to answer these fundamental questions by learning the continuous, time-dependent forces that guide cellular identity.

By framing cellular development as a Schrödinger Bridge problem, `scIDiff` goes beyond static snapshots of single-cell data. It reconstructs the complete "flow" of gene expression, revealing the dynamic regulatory landscape that cells navigate. This allows us to:

- **Learn the Regulatory Drift:** Uncover the vector field of gene regulation that pushes cells along developmental trajectories.
- **Quantify Irreversibility:** Measure the "point of no return" in cell fate decisions, identifying the moments when a cell commits to a specific lineage.
- **Infer Dynamic Gene Networks:** Move from static gene regulatory networks (GRNs) to a continuous movie of how gene-gene influences change over time.

`scIDiff` unifies the principles of optimal transport and score-based generative modeling to provide a physically-grounded, yet biologically-interpretable, model of cellular dynamics.

## 🌌 Mathematical Foundations

Given observed cell states sampled at different times or perturbation conditions, scIDiff models a stochastic process:

$$dX_t = u(X_t, t)dt + \sqrt{2\beta}dW_t$$

- $X_t \in \mathbb{R}^d$ — cell state (e.g., transcriptome)
- $u(X_t,t)$ — drift field (regulatory flow learned by scIDiff)
- $\beta$ — diffusion constant
- $dW_t$ — Brownian noise

The learned drift transports the empirical initial distribution $\rho_0$ to the terminal $\rho_1$ while minimizing a Schrödinger-Bridge energy:

$$\min_{u} \mathbb{E}[\int_0^1 \frac{1}{2\beta} \|u(X_t,t)\|^2 dt] \quad \text{s.t.} \quad X_0 \sim \rho_0, X_1 \sim \rho_1$$

The drift $u(x,t)$ satisfies the Schrödinger Bridge system:

$$u(x,t) = \nabla\log\phi(x,t) + \nabla\log\hat{\phi}(x,t)$$

where $(\phi,\hat{\phi})$ solve the Schrödinger half-bridge equations.

A temporal Jacobian

$$J(t) = \frac{\partial u}{\partial x}(t)$$

encodes local, causal gene-to-gene influence along time, where $J_{ij}(t) = \frac{\partial u_i}{\partial x_j}$ quantifies how gene $j$ influences the expression rate of gene $i$ at time $t$, providing dynamic GRN inference.

## 🔁 Forward and Reverse Drift

scIDiff learns both forward and reverse drift fields:

**Forward:**
$$dX_t = u(X_t, t)dt + \sqrt{2\beta}dW_t$$

**Reverse:**
$$dX_t = [u(X_t, t) - 2\beta\nabla_x\log\rho_t(X_t)]dt + \sqrt{2\beta}d\bar{W}_t$$

- The forward drift $u(x,t)$ describes how cells evolve naturally.
- The reverse drift captures how much "work" would be required to reprogram cells backward in time (e.g., iPSC or rejuvenation).

The difference $\Delta u = u_{\text{fwd}} - u_{\text{rev}}$ quantifies irreversibility — an analog of biological entropy production. This forward–reverse asymmetry provides a rigorous way to identify irreversible cell-fate decisions and the regulators (Yamanaka-like factors) capable of reversing them.

### Temporal Jacobians and Archetype Decomposition

scIDiff computes the temporal Jacobian tensors for both the forward ($J_{\text{fwd}}$) and reverse ($J_{\text{rev}}$) drift fields:

$$J_{\text{fwd}}(t) = \frac{\partial u_{\text{fwd}}}{\partial x}(t) \quad \text{and} \quad J_{\text{rev}}(t) = \frac{\partial u_{\text{rev}}}{\partial x}(t)$$

These Jacobians represent the instantaneous, causal influence of each gene on the expression rate of every other gene at a specific time $t$. To uncover the dominant modes of regulation and communication, scIDiff performs Singular Value Decomposition (SVD) on these tensors.

For a given Jacobian $J(t)$, the SVD is:

$$J(t) = U \Sigma V^T$$

- **Regulatory Archetypes (U):** The left singular vectors (columns of $U$) represent co-regulated gene modules, or *regulatory archetypes*, that are active during the process.
- **Communication Archetypes (V):** The right singular vectors (columns of $V$) correspond to the cellular states that are most influential, forming *communication archetypes*.

This decomposition provides a low-rank approximation of the complex, high-dimensional regulatory dynamics, revealing the underlying structure of gene regulation during both forward (e.g., differentiation) and reverse (e.g., reprogramming) processes.

### 🧬 Biological Guidance with RNA Velocity

While the Schrödinger Bridge finds the most probable path between two cell populations, this path is a mathematical optimum—the "straightest line" through a high-dimensional landscape. Biological processes, however, rarely take the straightest path. They are constrained by the intricate machinery of gene regulation.

To ensure our model respects these biological rules, `scIDiff` incorporates RNA velocity as a **biological compass**. By measuring the ratio of newly transcribed (unspliced) to mature (spliced) mRNA, RNA velocity provides a direct readout of the cell's immediate future—the direction it is *actually* moving in the next instant.

`scIDiff` uses this information not as a rigid constraint, but as a powerful guiding force. As it computes the long-range trajectory between cell states, it continuously checks its direction against the local, instantaneous vectors provided by RNA velocity. This ensures the global path is composed of steps that are mechanistically plausible and consistent with the cell's internal transcriptional dynamics.

The result is a model that balances global optimality with local biological realism, producing trajectories that are not just mathematically probable, but also biologically faithful.

## 🧠 What scIDiff Learns

| Layer                  | Description                                      | Output                            |
| ---------------------- | ------------------------------------------------ | --------------------------------- |
| Regulatory drift       | Deterministic direction of transcriptional change | $u(x,t)$                          |
| Reverse drift          | Reprogramming direction opposing natural evolution | $u_{\text{rev}}(x,t)$             |
| Irreversibility field  | Entropic asymmetry between forward and reverse   | $\Delta u = u - u_{\text{rev}}$   |
| Temporal Jacobian      | Local causal gene-to-gene influence              | $J(t)=\tfrac{\partial u}{\partial x}$ |
| Regulatory archetypes  | Low-rank modes of temporal regulation            | $u(x,t)\approx\sum_k a_k(t)A_kx$ |
| (Optional) Communication | Dynamic ligand–receptor signaling between cells    | $W_{ij}(t)$                       |

**Biological Insight:** The learned drift $u(x,t)$ captures:

- Attractors as zeros of $u(x,t)$ (fixed points)
- Differentiation paths as flow lines of $u(x,t)$
- Barriers as regions where $|u(x,t)|$ is large (rapid transitions)

## 🔗 Extension: scIDiff-Comm (Cell–Cell Communication)

With CellPhoneDB ligand–receptor priors, scIDiff models time-evolving communication graphs as coupled stochastic processes:

$$dX_t^{(i)} = u_{\text{intra}}(X_t^{(i)}, t)dt + \sum_j W_{ij}(t)f(X_t^{(j)} - X_t^{(i)})dt + \sqrt{2\beta}dW_t^{(i)}$$

- $W_{ij}(t)$ — communication strength between cells *i* and *j*
- $f(\Delta x)$ — interaction kernel (e.g., linear or gated)

**Result:** Communication archetypes (e.g., inflammatory relay, exhaustion/resolution) that co-evolve with regulatory drift.

## ⚖️ Quantifying Irreversibility

scIDiff computes entropy production and cycle flux metrics to measure the degree of biological irreversibility:

$$\dot{S}(t) = \frac{1}{\beta} \mathbb{E}[u(X_t,t) \cdot \nabla\log\rho_t(X_t)] = \mathbb{E}[\frac{\|u(X_t,t) - u_{\text{rev}}(X_t,t)\|^2}{2\beta}]$$

- High $\dot{S}(t)$ → irreversible differentiation (e.g., commitment, exhaustion)
- Low $\dot{S}(t)$ → reversible or plastic states (e.g., stem, progenitor)

This enables principled identification of control nodes capable of restoring reversibility.

## 🧮 Core Features

- Schrödinger-Bridge learning of forward & reverse dynamics
- Directional Jacobians for causal inference
- Quantitative irreversibility (entropy production, Δu)
- Counterfactual simulations for reprogramming paths
- Multi-omics ready (RNA, ATAC, velocity, protein, metabolite embeddings)
- Cytoscape exports for dynamic networks (.graphml, .cyjs)
- Optional scIDiff-Comm module for cell–cell signaling

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
```

## 🚀 Quick Start

```python
import scidiff
import scanpy as sc

# Load your single-cell data (AnnData format)
adata = sc.read_h5ad("your_data.h5ad")

# Initialize and train scIDiff
model = scidiff.SchrodingerBridge(adata, time_key="pseudotime")
model.train(n_epochs=1000)

# Extract regulatory dynamics
drift_field = model.get_drift()
jacobians = model.get_temporal_jacobians()
irreversibility = model.get_entropy_production()

# Simulate counterfactual reprogramming
reverse_paths = model.simulate_reverse_dynamics(n_cells=1000)
```

## 📚 Citation

If you use scIDiff in your research, please cite:

```bibtex
@article{scidiff2024,
  title={scIDiff: Schrödinger Bridge Learning of Single-Cell Regulatory Dynamics},
  author={Your Name and Collaborators},
  journal={Nature Methods},
  year={2024}
}
```

## 🤝 Contributing

We welcome contributions! Please see our Contributing Guidelines and Code of Conduct.

## 📄 License

This project is licensed under the MIT License - see the `LICENSE` file for details.

---

*scIDiff: Bridging cells across time through optimal transport and regulatory inference.*
