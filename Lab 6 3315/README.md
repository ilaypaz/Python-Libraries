Genomic Sequence Analysis Pipeline
Overview
This project implements a comprehensive bioinformatics pipeline for analyzing genomic sequences, focusing on sequence alignment, comparison, and visualization. The pipeline processes FASTA files containing nucleotide and amino acid sequences, performing tasks such as dot plot generation, pairwise and multiple sequence alignments, consensus sequence derivation, and phylogenetic analysis through dendrograms and tanglegrams. It is designed to support research in cancer genomics and other biological applications, aligning with clinical assay development goals.
Features

Data Retrieval: Automatically downloads FASTA files from specified URLs using urllib.request.
Sequence Analysis:
Computes sequence lengths for archaeal (Arch), selA, and HBA/HBB datasets.
Generates dot plots for sequence comparison with customizable window sizes and match thresholds using a custom dot_plot function.


Pairwise Alignments:
Performs global and local alignments for nucleotide and amino acid sequences using Biopython’s PairwiseAligner with BLOSUM62 and BLOSUM80 matrices.


Multiple Sequence Alignment (MSA):
Executes MSA using ClustalW for HBA/HBB nucleotide and amino acid sequences.
Derives consensus sequences for MSAs, handling ambiguous positions.


Visualization:
Creates publication-quality visualizations, including dot plots, MSA plots, dendrograms, and tanglegrams using Matplotlib and SciPy.
Supports segmented MSA visualization for specified position ranges.


Phylogenetic Analysis:
Computes distance matrices for aligned sequences.
Generates dendrograms and tanglegrams to compare sequence relationships across datasets.



Requirements

Python 3.x
Libraries: numpy, pandas, matplotlib, seaborn, biopython, scipy
External tool: ClustalW (clustalw2) installed and accessible in the system path
Operating System: Linux or compatible environment for file operations and ClustalW execution

Installation

Clone or download the project repository to your local machine.
Install required Python libraries:pip install numpy pandas matplotlib seaborn biopython scipy


Install ClustalW:
On Linux: sudo apt-get install clustalw (Ubuntu/Debian) or equivalent for your distribution.
Ensure clustalw2 is in your system path.


Set the working directory in the script to a valid path (e.g., replace /Users/ilaypaz/Downloads/Python with your directory).

Usage

Run the Script:
Execute the Python script (script.py) in a Linux environment:python script.py




Input Files:
The script downloads FASTA files from predefined URLs for archaeal (Arch), selA, and HBA/HBB sequences. Ensure internet connectivity or provide local FASTA files by modifying the URLs.


Output Files:
Dot Plots: Saved as PNG files (e.g., dotplot_w1_n1.png).
MSA Files: Saved as FASTA files (e.g., HBA_nuc_alignment.fasta).
Visualizations: Saved as PDF files (e.g., hba_aa_dendrogram.pdf, hba_aa_nuc_tanglegram.pdf).
Console Output: Includes sequence lengths, alignment scores, and consensus sequences.



Directory Structure

Input: FASTA files downloaded to the working directory (e.g., Arch.fa, selA.fa, HBA.fasta, HBA_aa.fasta).
Output:
Dot plots: PNG files in the working directory.
MSA alignments: FASTA files in the working directory.
Visualizations: PDF files for MSA plots, dendrograms, and tanglegrams.



Notes

The script assumes ClustalW is installed and accessible. If ClustalW fails, MSA-related steps (5b-6c) are skipped with a warning.
The pipeline processes HBA and HBB sequences identically in this implementation, as specified in the code.
Ensure sufficient disk space for temporary files created during ClustalW execution.
The script is optimized for a Linux environment; file path handling may need adjustment for other operating systems.

Limitations

Requires internet access for initial FASTA file downloads.
ClustalW dependency may introduce compatibility issues on non-Linux systems.
Visualization of large sequences may require significant memory and processing time.

Future Improvements

Add support for containerized execution (e.g., Docker) to manage dependencies like ClustalW.
Integrate Nextflow for workflow orchestration and scalability.
Enhance visualization with interactive outputs using Jupyter or Quarto notebooks.

Contact
For questions or support, contact the project supervisor or team lead, as this pipeline was developed for clinical bioinformatics applications similar to those at the Michael Smith Genome Sciences Centre.
