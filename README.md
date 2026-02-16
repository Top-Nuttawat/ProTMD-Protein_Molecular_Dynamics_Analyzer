# ProTMD - Protein Molecular Dynamics Analyzer

<img width="955" height="641" alt="image" src="https://github.com/user-attachments/assets/b7db94a4-8306-447e-b2ef-84afecffceda" />


**ProTMD** is a Python library designed for the comprehensive analysis of Molecular Dynamics (MD) trajectories. It bridges the gap between standard structural metrics and advanced graph-theoretical analysis, allowing researchers to explore protein flexibility, essential dynamics, and temporal residue-residue communication.

---

## 🚀 Key Features

- Structural Analytics: Automated computation of RMSD, RMSF, and Radius of Gyration ($R_g$).  

- Essential Dynamics (PCA): Perform Principal Component Analysis on $C\alpha$ trajectories to identify dominant collective motions.  

- Correlation Mapping: Dynamic Cross-Correlation Maps (DCCM) to visualize residue coupling.

- Temporal Protein Networks:
    - Dynamic Adjacency Matrices based on distance cutoffs.
    - Communicability Analysis: Uses matrix exponentials ($e^A$) to quantify all-atom communication pathways.

- Interactive Visualizations: High-quality plots using Plotly (3D surfaces, interactive line plots) and PyVis (interactive network graphs).

---

## 🛠 Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Top-Nuttawat/ProTMD-Protein_Molecular_Dynamics_Analyzer.git
cd ProTMD-Protein_Molecular_Dynamics_Analyzer
pip install -r requirements.txt
```
---
## 📖 Quick Start

1. Basic Structural Analysis


```bash
from ProTMD import Protein_Structure

# Initialize with GROMACS files
prot = Protein_Structure("protein.gro", "trajectory.xtc")

# Calculate and plot RMSD
prot.compute_RMSD()
prot.plot_RMSD(time_per_frame_ns=0.1)

# Principal Component Analysis
prot.compute_PCA()
prot.plot_PCA_scree(n_modes=10)
prot.write_eigenmode_XYZ(mode_index=0, amplitude=15.0) # Visualize PC1

```
2. Temporal Network Analysis

Explore how residues communicate over time using graph theory.

```bash
from ProTMD import Temporal_ProteinNetwork

# Build network with an 8.0Å cutoff
net = Temporal_ProteinNetwork("protein.gro", "trajectory.xtc", cutoff=8.0)

# Visualize the network at a specific frame (interactive HTML)
net.ProTNetworks(frame=100, html_name="residue_network.html")

# 3D Communicability Surface (exp(A))
net.plot3D_Communicability_matrix(frame=100)

```

---
## 📊 Module Overview

## `Protein_Structure`

Focuses on the geometry and collective motion of the protein.

| Method | Description |
|--------|------------|
| `compute_RMSD()` | Root Mean Square Deviation vs. Reference. |
| `compute_RMSF()` | Per-residue fluctuation (flexibility). |
| `compute_DCCM()` | Correlated motion between residue pairs. |
| `compute_PCA()` | Eigen-decomposition of the covariance matrix. |



##  `Temporal_ProteinNetwork`

Treats the protein as a dynamic graph where nodes are $C\alpha$ atoms.

| Method | Description |
|--------|------------|
| `Temporal_Adjacency_Matrix()` |Binary contact maps over time. |
| `Temporal_Communicability()` | Calculates $e^A$ to identify "hubs" of information flow.|
| `ProTNetworks()` | Exports interactive network graphs using NetworkX and PyVis. |

---

##  Dependencies
- MDAnalysis: Trajectory parsing
- NetworkX & PyVis: Graph theory and network visualization.
- Plotly: Interactive 3D and 2D plotting.
- Scikit-learn: PCA implementation.
- SciPy: Matrix exponentials for communicability.
---

##  ✍️ Author

Nuttawat Sawang  , M.Sc. (Physics)

Theoretical and Computational Physics (TCP) at Bangmod


📧 topza200915@gmail.com / nuttawatsawang.top@gmail.com

