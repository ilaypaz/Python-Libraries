HBA/HBB Phylogenetic Analysis
This repository contains a Python script for performing phylogenetic analysis on HBA (Hemoglobin Alpha) and HBB (Hemoglobin Beta) amino acid sequences. The script downloads FASTA files, conducts multiple sequence alignments (MSA) with ClustalW, constructs Neighbor-Joining (NJ) and Maximum Likelihood (ML) trees with bootstrap support, and generates cophylogenetic visualizations. This project demonstrates advanced bioinformatics techniques, automation, and data visualization skills.
Overview

Purpose: Analyze evolutionary relationships between HBA and HBB sequences across species using NJ and ML methods.
Features:
Downloads HBA and HBB amino acid sequences from a GitHub repository.
Performs MSA using ClustalW.
Constructs NJ and ML phylogenetic trees with 500+ bootstrap replicates.
Visualizes trees and cophylogenetic relationships using ete3 and matplotlib.
Converts scientific names to common names via NCBI Entrez API.
Compares clades and assesses taxonomic monophyly (e.g., snakes).


Date: Developed in Spring 2025.

Prerequisites

Python 3.8+
Required Python libraries:
biopython (for sequence handling and alignment)
numpy (for numerical operations)
pandas (for data manipulation)
matplotlib (for plotting)
ete3 (for tree manipulation and visualization)
requests (for HTTP requests)


External tools:
clustalw2 (for multiple sequence alignment)
phyml (for maximum likelihood tree construction, optional)


NCBI Entrez API:
Requires an email address for API access (set in the script).
Optional API key to avoid rate limits (obtain from NCBI).



Installation

Clone the Repository:
git clone https://github.com/ilaypaz/Bioinformatics.git
cd Bioinformatics


Set Up Virtual Environment (recommended):
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate


Install Dependencies:

Using pip:pip install biopython numpy pandas matplotlib ete3 requests


Using conda (recommended for bioinformatics tools):conda install -c bioconda clustalw phyml




Configure NCBI Entrez:

Edit the script (main.py) and replace your.email@example.com with your actual email.
Optional: Add an NCBI API key (commented in script) for faster access.


Verify Tools:

Ensure clustalw2 and phyml are in your PATH. Test with:clustalw2 -version
phyml --version





Usage

Run the Script:
python main.py


The script assumes the working directory is set to /Users/ilaypaz/Downloads/Python. Update os.chdir() in the script to your preferred directory.


Expected Output:

FASTA files: HBA_aa.fasta, HBB_aa.fasta
MSA files: HBA_aa_msa.fasta, HBB_aa_msa.fasta
Tree files: hba_aa_nj.tre, hbb_aa_nj.tre, hba_aa_ml.tre, hbb_aa_ml.tre
PDF plots: hba_aa_nj.pdf, hbb_aa_nj.pdf, hba_aa_ml.pdf, hbb_aa_ml.pdf, hba_aa_nj_bs.pdf, hbb_aa_nj_bs.pdf, hba_aa_ml_bs.pdf, hbb_aa_ml_bs.pdf, cophylo_njba_mlba.pdf, cophylo_njbb_mlbb.pdf, cophylo_mlba_mlbb.pdf
Console output: Sequence IDs, distance matrices, clade differences, and taxonomic analysis (e.g., snakes).


Notes:

The script includes error handling for missing dependencies (e.g., ete3, phyml).
NCBI rate limits may slow name conversion; use an API key if needed.
Temporary files are cleaned up automatically.



Outputs and Analysis

Phylogenetic Trees: NJ and ML trees are rooted at midpoint and visualized with bootstrap support, providing robust evolutionary insights.
Cophylogenetic Plots: Compare HBA vs. HBB and NJ vs. ML topologies, highlighting evolutionary congruence.
Clade Analysis: Identifies differences between HBA and HBB ML trees, and NJ vs. ML trees.
Taxonomic Monophyly: Checks for monophyletic groups (e.g., snakes) in HBA trees using NCBI taxonomy.

Acknowledgments

Data sourced from Biol-3315-files by Ido Hatam.
Utilizes BioPython, ete3, and NCBI Entrez API for bioinformatics analysis.
Inspired by coursework at Langara College, Spring 2025.

Contact

Author: Ilay Paz
Email: pazlilay@gmail.com
GitHub: github.com/ilaypaz

Last updated: May 19, 2025
