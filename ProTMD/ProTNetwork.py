import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

import networkx as nx
from pyvis.network import Network
from IPython.core.display import display, HTML

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

from scipy.linalg import expm


class Temporal_ProteinNetwork:
    def __init__(self , gro_file , xtc_file , cutoff =None , start_frame = None ):
        """
        Parameters:
        -----------
        gro_file : str
            Path to .gro file
        xtc_file : str
            Path to .xtc trajectory
        cutoff : float
            Distance cutoff (Å) for contact
        start_frame : int
            Frame to start network calculation
        """
        self.gro_file = gro_file
        self.xtc_file = xtc_file
        self.cutoff = cutoff
        self.start_frame = 0 if start_frame is None else start_frame

        # Load universe
        self.universe = mda.Universe(self.gro_file, self.xtc_file)
        self.CA_atoms = self.universe.select_atoms("name CA")
        self.n_residues = len(self.CA_atoms)
        self.n_frames = len(self.universe.trajectory) - self.start_frame

        self.CA_positions()
        self.Temporal_Adjacency_Matrix()
        self.Temporal_Communicability_Matrix()



    # Position of Cα atom for each residue over all trajectory frames
    def CA_positions(self):
        # Store coordinates for each frame
        CA_coords = []
        for ts in self.universe.trajectory:
            CA_coords.append(self.CA_atoms.positions.copy())

        self.CA_coords = np.array(CA_coords) # shape: (n_frames, n_CA, 3)
        print(f"Cα coordinates shape = {self.CA_coords.shape}")
        return self.CA_coords

    @staticmethod
    def distance( r1, r2 ):
        return ( (r1[0] - r2[0])**2 + (r1[1] - r2[1])**2 + (r1[2] - r2[2])**2 )**0.5

    def Temporal_Adjacency_Matrix(self):
        """Build temporal adjacency matrices A_time (frames x residues x residues)."""
        nstep = self.n_frames
        namino = self.n_residues
        self.A_time = np.zeros((nstep, namino, namino), dtype=int)
        self.pos = self.CA_coords  # shape: (nstep, namino, 3)
        print("Building temporal adjacency matrices")

        for istep in tqdm(range(nstep), desc="Frames"):
             for i in range(namino):
                for j in range(namino):
                    self.d = self.distance(self.pos[istep, i], self.pos[istep, j])
                    if self.d < self.cutoff and i != j:
                        self.A_time[istep, i, j] = 1

        print(f"Temporal adjacency matrix shape: {self.A_time.shape}")
        return self.A_time

    def ProTNetworks(self, frame=None, html_name="nx.html"):

        A = self.A_time[frame]
        G = nx.from_numpy_array(A)

        nt = Network('600px', '1000px', notebook=True, cdn_resources='in_line')
        nt.from_nx(G)
        nt.show(html_name)
        display(HTML(html_name))


    def plot_A_matrix(self, frame=None):
        plt.figure(figsize=(6, 6))
        plt.imshow(self.A_time[frame], interpolation='nearest')
        plt.title(f"Adjacency Matrix (frame {frame})")
        plt.colorbar()
        plt.show()

    def save_A_time(self, filename=None):
        base = os.path.splitext(os.path.basename(self.gro_file))[0]

        if filename is None:
            filename = f"{base}-Atime-r{self.cutoff}.dat"
        
        # Flatten to 1D
        A_time_1d = self.A_time.reshape(-1)
        np.savetxt(filename, A_time_1d, fmt='%1d')
        print(f"Saved A_time → {filename}")

    #The Communicability Matrix

    def  Temporal_Communicability_Matrix(self):
        """
        Compute G_matrix_function_time[t] = expm(A_time[t])

        """
        n_frames, namino, _ = self.A_time.shape
        self.G_matrix_function_time = np.zeros((n_frames, namino, namino))

        for i in tqdm(range(n_frames), desc="Computing expm(A_time[t])"):
            self.G_matrix_function_time[i] = expm(self.A_time[i, :, :])

        print(f"G_matrix_function_time shape: {self.G_matrix_function_time.shape}")
        return self.G_matrix_function_time


    def plotmap_Communicability_Matrix(self, frame=None):
        plt.figure(figsize=(6, 6))
        plt.imshow(self.G_matrix_function_time[frame], interpolation='nearest')
        plt.title(f"Communicability_Matrix (frame {frame})")
        plt.colorbar()
        plt.show()

    def plot3D_Communicability_matrix(self, frame=None):
        """Plot exp(A) as a 3D surface for selected frame."""

        self.Temporal_Communicability_Matrix()

        A_exp = self.G_matrix_function_time[frame]

        fig = go.Figure(
            data=[go.Surface(z=A_exp)]
        )

        fig.update_layout(
            title=f'Matrix Exponential 3D Surface Map (frame {frame})',
            scene=dict(
                xaxis_title='Node j',
                yaxis_title='Node i',
                zaxis_title='exp(A)'
            ),
            width=1000,
            height=1000
        )

        fig.show()


    def save_G_time(self, filename=None):
        base = os.path.splitext(os.path.basename(self.gro_file))[0]

        if filename is None:
            filename = f"{base}-Gtime-r{self.cutoff}.dat"

        # Flatten to 1D
        G_time_1d = self.G_matrix_function_time.reshape(-1)

        # Save as float
        np.savetxt(filename, G_time_1d, fmt='%.6f')

        print(f"Saved G_time → {filename}")
