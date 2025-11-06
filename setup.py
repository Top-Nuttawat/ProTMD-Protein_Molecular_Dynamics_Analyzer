from setuptools import setup, find_packages

setup(
    name="ProTMD",
    version="0.1.0",
    author="Nuttawat Sawang",
    author_email="topza200915@gmail.com",
    description="Protein Trajectory Analysis and Essential Dynamics tools using MDAnalysis",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Top-Nuttawat/ProTMD-Protein_Molecular_Dynamics_Analyzer",
    packages=find_packages(),
    install_requires=open("requirements.txt").read().splitlines(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Molecular BioPhysics and Bio Informatics",
    ],
    python_requires=">=3.8",
)
