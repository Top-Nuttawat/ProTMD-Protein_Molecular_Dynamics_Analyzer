import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

import os

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

from sklearn.decomposition import PCA


class Protein_Structure:
    
    """
    ProTMD: Protein Trajectory Motion Dynamics Library
    Provides RMSD, RMSF, Rg, DCCM, and PCA analysis from MD trajectories.
    
    """
    
    def __init__(self, gro_file, xtc_file):
        self.gro_file = gro_file
        self.xtc_file = xtc_file
        self.universe = mda.Universe(self.gro_file, self.xtc_file)
        self.CA_positions()


# Position of Cα atom for each residue over all trajectory frames
    def CA_positions(self):
       # Select all alpha carbon atoms for each frame
        self.CA_atoms = self.universe.select_atoms("name CA")

        # Store coordinates for each frame
        CA_coords = []
        for ts in self.universe.trajectory:
            CA_coords.append(self.CA_atoms.positions.copy())

        self.CA_coords = np.array(CA_coords) # shape: (n_frames, n_CA, 3)
        print(f"Cα coordinates shape = {self.CA_coords.shape}")
        return self.CA_coords


    #Write the Cα coordinates to an XYZ file
    def write_CA_XYZ(self, output_filename=None):
        if output_filename is None:
            if '.' in self.gro_file:
                base_name = self.gro_file.split('.')[0]
            else:
                base_name = self.gro_file
            output_filename = base_name + "_c_alpha_trajectory.xyz"

        with open(output_filename, 'w') as outfile:
            for timestep_coords in self.CA_coords:
                num_atoms = timestep_coords.shape[0]
                outfile.write(f"{num_atoms}\n")
                outfile.write("Generated from MDAnalysis Cα trajectory\n")
                for atom_coord in timestep_coords:
                    outfile.write(f"CA {atom_coord[0]:.3f} {atom_coord[1]:.3f} {atom_coord[2]:.3f}\n")

        print(f"Cα trajectory written to {output_filename}")


    # RMSD
    def compute_RMSD(self, reference_frame=0):
        ref = self.CA_coords[reference_frame]
        num_atoms = ref.shape[0]
        self.rmsd = np.sqrt(np.mean(np.sum((self.CA_coords - ref)**2, axis=2), axis=1))
        return self.rmsd

    def plot_RMSD(self, time_per_frame_ns=0.1):
        time_ns = np.arange(len(self.rmsd)) * time_per_frame_ns
        fig = px.line(
            x=time_ns,
            y=self.rmsd,
            labels={"x": "Time (ns)", "y": "RMSD (Å)"},
            title="Cα RMSD "
        )
        fig.update_traces(line_color="blue")
        fig.update_layout(template="plotly_white")
        fig.show()

    # RMSF
    def compute_RMSF(self):
        mean_coords = np.mean(self.CA_coords, axis=0)
        diff = self.CA_coords - mean_coords
        self.rmsf = np.sqrt(np.mean(np.sum(diff**2, axis=2), axis=0))
        return self.rmsf

    def plot_RMSF(self):
        fig = px.line(
            x=np.arange(len(self.rmsf)) + 1,
            y=self.rmsf,
            labels={"x": "Residue index", "y": "RMSF (Å)"},
            title="Cα RMSF per Residue"
        )
        fig.update_traces(line_color="green")
        fig.update_layout(template="plotly_white")
        fig.show()

    # Radius of gyration
    def compute_Rg(self):
        # Center of geometry per frame
        com = self.CA_coords.mean(axis=1, keepdims=True)

        # Radius of gyration per frame
        self.Rg = np.sqrt(np.mean(np.sum((self.CA_coords - com)**2, axis=2), axis=1))
        return self.Rg

    def plot_Rg(self, time_per_frame_ns=0.1):
        time_ns = np.arange(len(self.Rg)) * time_per_frame_ns
        fig = px.line(
            x=time_ns,
            y=self.Rg,
            labels={"x": "Time (ns)", "y": "Rg (Å)"},
            title="Radius of gyration over Time"
        )
        fig.update_traces(line_color="purple")
        fig.update_layout(template="plotly_white")
        fig.show()

    # DCCM
    def compute_DCCM(self):
        """
        Compute Dynamic Cross-Correlation Map (DCCM)
        """
        # Mean-subtracted coordinates per residue
        centered_coords = self.CA_coords - np.mean(self.CA_coords, axis=0)  # (n_frames, n_CA, 3)
        n_frames, n_res, _ = centered_coords.shape

        # Initialize covariance matrix
        covariance_matrix = np.zeros((n_res, n_res))

        # Compute covariance
        for i in range(n_res):
            for j in range(n_res):
                dot_products = np.sum(centered_coords[:, i, :] * centered_coords[:, j, :], axis=1)
                covariance_matrix[i, j] = np.mean(dot_products)

        # Standard deviation per residue
        std_devs = np.sqrt(np.diag(covariance_matrix))
        std_devs[std_devs == 0] = 1e-10  # avoid division by zero

        # Normalize to get DCCM
        self.DCCM = covariance_matrix / np.outer(std_devs, std_devs)
        return self.DCCM


    def plot_DCCM(self):
        """
        Plot DCCM using Plotly heatmap.
        """
        fig = px.imshow(
            self.DCCM,
            labels=dict(x="Residue", y="Residue", color="Correlation"),
            x=np.arange(1, self.DCCM.shape[0]+1),
            y=np.arange(1, self.DCCM.shape[0]+1),
            color_continuous_scale='RdBu',
            zmin=-1, zmax=1,
            title="Dynamic Cross-Correlation Map (DCCM)"
        )
        fig.update_layout(template="plotly_white")
        fig.show()

    # Principal Component Analysis (Essential Dynamics)
    def compute_PCA(self, start=0, end=None, trim_residues=None):
        """
        Perform PCA (Essential Dynamics) on Cα trajectory.
        Parameters:
        -----------
        start : int
            Starting frame for PCA analysis
        end : int
            Ending frame for PCA analysis
        trim_residues : tuple or list (optional)
            Range of residues to include (e.g., (0, 375) to exclude H8 loop)
        """
        # 1. Select frames
        coords = self.CA_coords[start:end]  # shape: (n_frames, n_res, 3)

        # 2. Optionally trim residues (e.g., remove H8 loop)
        if trim_residues is not None:
            coords = coords[:, trim_residues[0]:trim_residues[1], :]

        # 3. Compute average structure
        avg_coords = np.mean(coords, axis=0)

        # 4. Center the trajectory
        centered = coords - avg_coords

        # 5. Reshape into 2D (frames × (res*3))
        X = centered.reshape(centered.shape[0], -1)

        # 6. Compute covariance matrix
        cov = np.cov(X.T)

        # 7. Eigen decomposition
        eigenvalues, eigenvectors = np.linalg.eigh(cov)

        # 8. Sort eigenmodes by descending eigenvalue
        sort_idx = np.argsort(eigenvalues)[::-1]
        self.eigenvalues = eigenvalues[sort_idx]
        self.eigenvectors = eigenvectors[:, sort_idx]
        self.avg_coords = avg_coords
        self.trim_shape = coords.shape[1:]  # (n_res, 3)

        print(" PCA complete.")
        print(f"Shape of eigenvectors: {self.eigenvectors.shape}")
        print(f"Top 5 eigenvalues: {self.eigenvalues[:5]}")

        return self.eigenvalues, self.eigenvectors


    def plot_PCA_scree(self, n_modes=20):
        """
        Plot eigenvalue spectrum (scree plot)
        """
        total_var = np.sum(self.eigenvalues)
        var_exp = (self.eigenvalues / total_var) * 100
        cum_var = np.cumsum(var_exp)

        fig = px.bar(
            x=np.arange(1, n_modes+1),
            y=var_exp[:n_modes],
            labels={"x": "Principal Component", "y": "Variance explained (%)"},
            title="PCA Scree Plot (Essential Dynamics)"
        )
        fig.add_scatter(x=np.arange(1, n_modes+1), y=cum_var[:n_modes], mode="lines+markers", name="Cumulative")
        fig.update_layout(template="plotly_white")
        fig.show()


    def write_eigenmode_XYZ(self, mode_index=None , amplitude=10.0, num_frames=50, output_filename=None):
        """
        Generate and write an XYZ trajectory visualizing motion along a given eigenmode.
        - Choose the eigenmode to visualize (e.g., the first mode, index 0)
        - Choose an amplitude for scaling the motion along the eigenmode
        - A larger amplitude will exaggerate the motion for visualization
        -  Generate frames along the eigenmode num_frames = 50  Number of frames to generate for the trajectory
        """
        # Choose eigenvector
        eigvec = self.eigenvectors[:, mode_index]
        eigvec_reshaped = eigvec.reshape(self.trim_shape)  # Reshape the eigenvector to have dimensions (num_atoms, 3)

        frames = []
        for i in range(num_frames):
            displacement = amplitude * np.sin(2 * np.pi * i / (num_frames - 1)) * eigvec_reshaped
            frame_coords = self.avg_coords + displacement
            frames.append(frame_coords)
        eigenmode_trajectory = np.array(frames)

        # Output filename
        if output_filename is None:
            base_name = os.path.splitext(os.path.basename(self.gro_file))[0]
            output_filename = f"{base_name }_eigenmode_{mode_index+1}_trajectory.xyz"

        with open(output_filename, "w") as outfile:
            for timestep_coords in eigenmode_trajectory:
                num_atoms = timestep_coords.shape[0]
                outfile.write(f"{num_atoms}\n")
                outfile.write(f"Eigenmode {mode_index+1} (amplitude {amplitude})\n")
                for atom_coord in timestep_coords:
                    outfile.write(f"CA {atom_coord[0]:.3f} {atom_coord[1]:.3f} {atom_coord[2]:.3f}\n")

        print(f" Eigenmode trajectory written to {output_filename}")
        return output_filename , eigvec_reshaped