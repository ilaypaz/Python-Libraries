import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Align.Applications import ClustalwCommandline
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import os
import tempfile
import re
from io import StringIO

# Set working directory
os.chdir("/Users/ilaypaz/Downloads/Python")

# Custom dot plot function (emulating seqinr::dotPlot)
def dot_plot(seq1, seq2, wsize=1, nmatch=1, col=("black", "white")):
    """
    Create a dot plot comparing two sequences.
    Parameters:
    - seq1, seq2: Sequences (str or Bio.Seq)
    - wsize: Window size for comparison
    - nmatch: Number of matches required in window
    - col: Colors for matches (match, no-match)
    """
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    len1, len2 = len(seq1), len(seq2)
    matrix = np.zeros((len1, len2))
    
    for i in range(len1 - wsize + 1):
        for j in range(len2 - wsize + 1):
            window1 = seq1[i:i+wsize]
            window2 = seq2[j:j+wsize]
            matches = sum(a == b for a, b in zip(window1, window2))
            if matches >= nmatch:
                matrix[i, j] = 1
    
    plt.figure(figsize=(8, 8))
    plt.imshow(matrix, cmap=plt.cm.colors.ListedColormap(col), interpolation='nearest')
    plt.xlabel("Sequence 2")
    plt.ylabel("Sequence 1")
    plt.title(f"Dot Plot: wsize={wsize}, nmatch={nmatch}")
    plt.savefig(f"dotplot_w{wsize}_n{nmatch}.png")
    plt.close()

# Download and read FASTA files
def download_fasta(url, filename):
    try:
        urllib.request.urlretrieve(url, filename)
        return list(SeqIO.parse(filename, "fasta"))
    except Exception as e:
        print(f"Error downloading {filename}: {e}")
        return []

arch_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/Arch.fa"
selA_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/selA.fa"
hba_nuc_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA.fasta"
hba_aa_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta"

arch = download_fasta(arch_url, "Arch.fa")
selA = download_fasta(selA_url, "selA.fa")
nucs = download_fasta(hba_nuc_url, "HBA.fasta")
aas = download_fasta(hba_aa_url, "HBA_aa.fasta")

if not all([arch, selA, nucs, aas]):
    print("Error: One or more FASTA files failed to load. Check URLs or download manually.")
    exit(1)

# 1a: Sequence lengths and dot plots for arch
print("Arch sequence lengths:", [len(seq.seq) for seq in arch])
dot_plot(arch[0].seq, arch[2].seq, col=("black", "white"))
dot_plot(arch[1].seq, arch[2].seq, col=("black", "white"))

# 1c: Dot plots with window size and matches
dot_plot(arch[0].seq, arch[2].seq, wsize=50, nmatch=30)
dot_plot(arch[1].seq, arch[2].seq, wsize=50, nmatch=30)

# 2a: selA sequence lengths and dot plot
print("selA sequence lengths:", [len(seq.seq) for seq in selA])
dot_plot(selA[0].seq, selA[4].seq, col=("black", "white"))

# 2b: Dot plot with wsize=10, nmatch=5
dot_plot(selA[0].seq, selA[4].seq, wsize=10, nmatch=5)

# 2c: Dot plot with wsize=60, nmatch=12
dot_plot(selA[0].seq, selA[4].seq, wsize=60, nmatch=12)

# 3a: Pairwise alignments
# Prepare arch_16s (last 3 sequences)
arch_16s = [str(seq.seq) for seq in arch[-3:]]
selAset = [str(seq.seq).upper() for seq in selA]

# Nucleotide alignments
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1

globalseq13 = aligner.align(arch_16s[0], arch_16s[2])[0]
globalseq23 = aligner.align(arch_16s[1], arch_16s[2])[0]
print("Global seq 1-3 score:", globalseq13.score)
print("Global seq 2-3 score:", globalseq23.score)

# Amino acid alignments
aligner_aa = PairwiseAligner()
aligner_aa.mode = "global"

# Available matrices
available_matrices = substitution_matrices.load()
print("Available substitution matrices:", available_matrices)

# Dictionary to store alignments
aa_alignments = {}

# BLOSUM62 (default)
try:
    aligner_aa.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aa_alignments["BLOSUM62"] = aligner_aa.align(selAset[1], selAset[2])[0]
except Exception as e:
    print(f"Error loading BLOSUM62: {e}")

# BLOSUM80
try:
    aligner_aa.substitution_matrix = substitution_matrices.load("BLOSUM80")
    aa_alignments["BLOSUM80"] = aligner_aa.align(selAset[1], selAset[2])[0]
except Exception as e:
    print(f"Error loading BLOSUM80: {e}")

print("AA alignments (selA 2-3):")
for matrix, alignment in aa_alignments.items():
    print(f"{matrix} score:", alignment.score)

# 4a: Local alignments
aligner.mode = "local"
localseq13 = aligner.align(arch_16s[0], arch_16s[2])[0]
localseq23 = aligner.align(arch_16s[1], arch_16s[2])[0]

aligner_aa.mode = "local"
local_aa_alignments = {}

try:
    aligner_aa.substitution_matrix = substitution_matrices.load("BLOSUM62")
    local_aa_alignments["BLOSUM62"] = aligner_aa.align(selAset[0], selAset[4])[0]
except Exception as e:
    print(f"Error loading BLOSUM62 for local: {e}")

try:
    aligner_aa.substitution_matrix = substitution_matrices.load("BLOSUM80")
    local_aa_alignments["BLOSUM80"] = aligner_aa.align(selAset[0], selAset[4])[0]
except Exception as e:
    print(f"Error loading BLOSUM80 for local: {e}")

# 4b: Analyze local alignment (seq 1-3)
aligned_pattern = str(localseq13.aligned[0])
aligned_subject = str(localseq13.aligned[1])
aligned_pattern_no_gaps = aligned_pattern.replace("-", "")
aligned_subject_no_gaps = aligned_subject.replace("-", "")
print("Local seq 1-3 pattern length (no gaps):", len(aligned_pattern_no_gaps))
print("Local seq 1-3 subject length (no gaps):", len(aligned_subject_no_gaps))

# 4c: Verify archaea (manual check in R; here we print description)
print("Arch[1] description:", arch[0].description)
# Note: Manual verification confirms Aeropyrum camini (archaea)

# 4d: AA alignment lengths
print("AA alignment lengths (selA 2-3):")
for matrix, alignment in aa_alignments.items():
    print(f"{matrix} length:", len(str(alignment.aligned[0]).replace("-", "")))

# 5a: Process HBA/HBB names
# Nucleotide sequences
nuc_names = [re.split(r"\s+", rec.description)[1:3] for rec in nucs]
nuc_names = ["_".join(name).replace(" ", "_") for name in nuc_names]
# Ensure unique names by appending index if needed
seen = {}
for i, name in enumerate(nuc_names):
    if name in seen:
        seen[name] += 1
        nuc_names[i] = f"{name}_{seen[name]}"
    else:
        seen[name] = 1

# Amino acid sequences
aa_names = [re.split(r"\[|\]", rec.description)[1] for rec in aas]
aa_names = [name.replace(" ", "_") for name in aa_names]
# Ensure unique names by appending index
seen = {}
for i, name in enumerate(aa_names):
    if name in seen:
        seen[name] += 1
        aa_names[i] = f"{name}_{seen[name]}"
    else:
        seen[name] = 1

# Assign names and print for debugging
print("Nucleotide sequence IDs:", nuc_names)
print("Amino acid sequence IDs:", aa_names)

for i, rec in enumerate(nucs):
    rec.id = nuc_names[i]
    rec.description = ""
for i, rec in enumerate(aas):
    rec.id = aa_names[i]
    rec.description = ""

# 5b: MSA with ClustalW
def run_clustalw(input_seqs, output_file):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
            SeqIO.write(input_seqs, tmp.name, "fasta")
            clustalw_cline = ClustalwCommandline("clustalw2", infile=tmp.name, outfile=output_file, output="FASTA")
            stdout, stderr = clustalw_cline()
        return AlignIO.read(output_file, "fasta")
    except Exception as e:
        print(f"Error running ClustalW: {e}")
        return None

# Initialize alignments
hbanucs = run_clustalw(nucs, "HBA_nuc_alignment.fasta")
hbaaas = run_clustalw(aas, "HBA_aa_alignment.fasta")
hbbnucs = run_clustalw(nucs, "HBB_nuc_alignment.fasta")  # Same as HBA for now
hbbaas = run_clustalw(aas, "HBB_aa_alignment.fasta")

# Check if alignments succeeded
if not all([hbanucs, hbaaas, hbbnucs, hbbaas]):
    print("Warning: One or more ClustalW alignments failed. Skipping MSA-related steps (5b-6c).")
else:
    # 5c: Consensus sequences
    def msa_consensus(alignment):
        seqs = [str(rec.seq) for rec in alignment]
        length = len(seqs[0])
        consensus = ""
        for i in range(length):
            column = [seq[i] for seq in seqs]
            most_common = max(set(column), key=column.count)
            consensus += most_common if column.count(most_common) > len(seqs) / 2 else "?"
        return consensus.replace("?", "-") if "nuc" in str(alignment[0].seq).lower() else consensus.replace("?", "*")

    hba_n_msa_consensus = msa_consensus(hbanucs)
    hba_aa_msa_consensus = msa_consensus(hbaaas)
    hbb_n_msa_consensus = msa_consensus(hbbnucs)
    hbb_aa_msa_consensus = msa_consensus(hbbaas)

    # Append consensus to alignments
    hba_n_msa_ms = hbanucs[:]
    hba_n_msa_ms.append(SeqRecord(Seq(hba_n_msa_consensus), id="Consensus"))
    hba_aa_msa_ms = hbaaas[:]
    hba_aa_msa_ms.append(SeqRecord(Seq(hba_aa_msa_consensus), id="Consensus"))
    hbb_n_msa_ms = hbbnucs[:]
    hbb_n_msa_ms.append(SeqRecord(Seq(hbb_n_msa_consensus), id="Consensus"))
    hbb_aa_msa_ms = hbbaas[:]
    hbb_aa_msa_ms.append(SeqRecord(Seq(hbb_aa_msa_consensus), id="Consensus"))

    SeqIO.write(hba_n_msa_ms, "HBA_nuc_alignment.fasta", "fasta")
    SeqIO.write(hba_aa_msa_ms, "HBA_aa_alignment.fasta", "fasta")
    SeqIO.write(hbb_n_msa_ms, "HBB_nuc_alignment.fasta", "fasta")
    SeqIO.write(hbb_aa_msa_ms, "HBB_aa_alignment.fasta", "fasta")

    print("HBA nuc consensus:", hba_n_msa_consensus)
    print("HBA aa consensus:", hba_aa_msa_consensus)
    print("HBB nuc consensus:", hbb_n_msa_consensus)
    print("HBB aa consensus:", hbb_aa_msa_consensus)

    # 5d: MSA visualizations (approximating ggmsa)
    def plot_msa(alignment, start, end, filename, title):
        seqs = [str(rec.seq)[start-1:end] for rec in alignment]
        plt.figure(figsize=(10, len(seqs) * 0.5))
        for i, seq in enumerate(seqs):
            for j, char in enumerate(seq):
                plt.text(j + 0.5, i + 0.5, char, ha="center", va="center", fontsize=8)
        plt.xlim(0, end - start + 1)
        plt.ylim(0, len(seqs))
        plt.title(title)
        plt.xlabel("Position")
        plt.ylabel("Sequence")
        plt.savefig(filename)
        plt.close()

    for i in range(1, 144, 50):
        end = min(i + 49, len(hba_aa_msa_ms[0].seq))
        plot_msa(hba_aa_msa_ms, i, end, f"HBA_alignment_positions_{i}_to_{end}.pdf",
                 f"HBA Alignment Positions {i} to {end}")

    for i in range(1, 148, 50):
        end = min(i + 49, len(hbb_aa_msa_ms[0].seq))
        plot_msa(hbb_aa_msa_ms, i, end, f"HBB_alignment_positions_{i}_to_{end}.pdf",
                 f"HBB Alignment Positions {i} to {end}")

    # 5e: Distance matrices
    def dist_alignment(alignment):
        seqs = [str(rec.seq) for rec in alignment]
        n = len(seqs)
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                pairs = sum(a != b and a != "-" and b != "-" for a, b in zip(seqs[i], seqs[j]))
                dist[i, j] = dist[j, i] = pairs / len(seqs[i].replace("-", ""))
        return squareform(dist)  # Convert to condensed format

    hba_aa_msa_m_dist = dist_alignment(hba_aa_msa_ms)
    hba_n_msa_m_dist = dist_alignment(hba_n_msa_ms)
    hbb_aa_msa_m_dist = dist_alignment(hbb_aa_msa_ms)
    hbb_n_msa_m_dist = dist_alignment(hbb_n_msa_ms)

    # 5f: Dendrograms
    def plot_dendrogram(dist, title, filename):
        linkage_matrix = linkage(dist, method="average")
        plt.figure(figsize=(8, 6))
        dendrogram(linkage_matrix, labels=[rec.id for rec in hba_aa_msa_ms])
        plt.title(title)
        plt.savefig(filename)
        plt.close()

    plot_dendrogram(hba_aa_msa_m_dist, "HBA AA Dendrogram", "hba_aa_dendrogram.pdf")
    plot_dendrogram(hba_n_msa_m_dist, "HBA Nuc Dendrogram", "hba_nuc_dendrogram.pdf")
    plot_dendrogram(hbb_aa_msa_m_dist, "HBB AA Dendrogram", "hbb_aa_dendrogram.pdf")
    plot_dendrogram(hbb_n_msa_m_dist, "HBB Nuc Dendrogram", "hbb_nuc_dendrogram.pdf")

    # 6a-c: Tanglegrams (approximating dendextend::tanglegram)
    def plot_tanglegram(dist1, dist2, title, filename, labels1, labels2):
        linkage1 = linkage(dist1, method="average")
        linkage2 = linkage(dist2, method="average")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        dend1 = dendrogram(linkage1, ax=ax1, labels=labels1)
        dend2 = dendrogram(linkage2, ax=ax2, labels=labels2, orientation="right")
        
        # Connect matching labels
        for i, label1 in enumerate(dend1["ivl"]):
            for j, label2 in enumerate(dend2["ivl"]):
                if label1 == label2:
                    ax1.plot([1, 1.5], [i, i], "k-", alpha=0.3)
                    ax2.plot([-1, -1.5], [j, j], "k-", alpha=0.3)
                    plt.plot([ax1.get_xlim()[1], ax2.get_xlim()[0]], [i, j], "k-", alpha=0.3)
        
        ax1.set_title(f"{title} (Left: First Dendrogram)")
        ax2.set_title(f"{title} (Right: Second Dendrogram)")
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()

    plot_tanglegram(hba_aa_msa_m_dist, hba_n_msa_m_dist, "HBA AA vs Nuc Tanglegram",
                    "hba_aa_nuc_tanglegram.pdf", [rec.id for rec in hba_aa_msa_ms], [rec.id for rec in hba_n_msa_ms])
    plot_tanglegram(hbb_aa_msa_m_dist, hbb_n_msa_m_dist, "HBB AA vs Nuc Tanglegram",
                    "hbb_aa_nuc_tanglegram.pdf", [rec.id for rec in hbb_aa_msa_ms], [rec.id for rec in hbb_n_msa_ms])
    plot_tanglegram(hbb_aa_msa_m_dist, hba_aa_msa_m_dist, "HBB AA vs HBA AA Tanglegram",
                    "hbb_aa_hba_aa_tanglegram.pdf", [rec.id for rec in hbb_aa_msa_ms], [rec.id for rec in hba_aa_msa_ms])
