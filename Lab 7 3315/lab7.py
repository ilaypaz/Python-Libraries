import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
import re
import time
try:
    from ete3 import Tree, TreeStyle, TextFace
except ImportError:
    print("Error: ete3 is not installed. Install with: pip install ete3")
    exit(1)
import requests
from io import StringIO

# Set working directory
os.chdir("/Users/ilaypaz/Downloads/Python")

# NCBI Entrez email (required for API access)
Entrez.email = "your.email@example.com"  # REPLACE WITH YOUR EMAIL (e.g., ilaypaz@example.com)
# Optional: Add NCBI API key to avoid rate limits
# Entrez.api_key = "your_api_key"  # Obtain from NCBI

# Function to download FASTA files
def download_fasta(url, filename):
    try:
        urllib.request.urlretrieve(url, filename)
        return list(SeqIO.parse(filename, "fasta"))
    except Exception as e:
        print(f"Error downloading {filename}: {e}")
        return []

# Load HBA and HBB amino acid sequences
hba_aa_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta"
hbb_aa_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBB_aa.fasta"

hba_aa = download_fasta(hba_aa_url, "HBA_aa.fasta")
hbb_aa = download_fasta(hbb_aa_url, "HBB_aa.fasta")

if not (hba_aa and hbb_aa):
    print("Error: Failed to load FASTA files. Exiting.")
    exit(1)

# Clean names (mimic R's str_split_fixed and str_remove_all)
def clean_names(records):
    names = []
    for rec in records:
        match = re.search(r"\[(.*?)\]", rec.description)
        name = match.group(1) if match else rec.id
        names.append(name.replace(" ", "_"))
    return names

hba_names = clean_names(hba_aa)
hbb_names = clean_names(hbb_aa)

# Ensure both datasets have the same organisms
common_names = set(hba_names).intersection(hbb_names)
hba_aa = [rec for rec, name in zip(hba_aa, hba_names) if name in common_names]
hbb_aa = [rec for rec, name in zip(hbb_aa, hbb_names) if name in common_names]

# Update names in records
for i, (rec, name) in enumerate(zip(hba_aa, [n for n in hba_names if n in common_names])):
    rec.id = name
    rec.description = ""
for i, (rec, name) in enumerate(zip(hbb_aa, [n for n in hbb_names if n in common_names])):
    rec.id = name
    rec.description = ""

print("HBA sequence IDs:", [rec.id for rec in hba_aa])
print("HBB sequence IDs:", [rec.id for rec in hbb_aa])

# Perform MSA with ClustalW
def run_clustalw(input_seqs, output_file):
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
            SeqIO.write(input_seqs, tmp.name, "fasta")
            clustalw_cline = ClustalwCommandline("clustalw2", infile=tmp.name, outfile=output_file, output="FASTA")
            stdout, stderr = clustalw_cline()
        return AlignIO.read(output_file, "fasta")
    except Exception as e:
        print(f"Error running ClustalW: {e}. Ensure ClustalW is installed (conda install -c bioconda clustalw)")
        return None

hba_aa_msa = run_clustalw(hba_aa, "HBA_aa_msa.fasta")
hbb_aa_msa = run_clustalw(hbb_aa, "HBB_aa_msa.fasta")

if not (hba_aa_msa and hbb_aa_msa):
    print("Error: MSA failed. Exiting.")
    exit(1)

# Check if HBA and HBB sequences are identical
print("HBA AA seqs equal HBB AA seqs:", all(str(a.seq) == str(b.seq) for a, b in zip(hba_aa_msa, hbb_aa_msa)))

# Generate distance matrices (approximating JTT model)
calculator = DistanceCalculator('identity')
hba_aa_dml = calculator.get_distance(hba_aa_msa)
hbb_aa_dml = calculator.get_distance(hbb_aa_msa)

print("HBA distance matrix:\n", hba_aa_dml.matrix)
print("HBB distance matrix:\n", hbb_aa_dml.matrix)

# Generate NJ trees
try:
    constructor = DistanceTreeConstructor(calculator)
    hba_aa_nj = constructor.nj(hba_aa_dml)
    hbb_aa_nj = constructor.nj(hbb_aa_dml)
except Exception as e:
    print(f"Error building NJ trees: {e}")
    exit(1)

# Root trees at midpoint
def midpoint_root(tree):
    try:
        newick_str = tree.format('newick')
        ete_tree = Tree(newick_str, format=1)  # Use format=1 to handle Biopython's internal node labels
        outgroup = ete_tree.get_midpoint_outgroup()
        if outgroup:
            ete_tree.set_outgroup(outgroup)
        return ete_tree
    except Exception as e:
        print(f"Error rooting tree: {e}. Newick string: {newick_str}")
        return None

hba_aa_nj_ete = midpoint_root(hba_aa_nj)
hbb_aa_nj_ete = midpoint_root(hbb_aa_nj)

if not (hba_aa_nj_ete and hbb_aa_nj_ete):
    print("Error: Failed to root NJ trees. Exiting.")
    exit(1)

# Save NJ trees
hba_aa_nj_ete.write(outfile="hba_aa_nj.tre", format=0)
hbb_aa_nj_ete.write(outfile="hbb_aa_nj.tre", format=0)

# Plot NJ trees
def plot_tree(tree, title, filename):
    try:
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.title.add_face(TextFace(title, fsize=10), column=0)
        tree.render(filename, tree_style=ts)
    except Exception as e:
        print(f"Error plotting tree {filename}: {e}")

plot_tree(hba_aa_nj_ete, "Rooted Neighbor Joining Tree of HBA", "hba_aa_nj.pdf")
plot_tree(hbb_aa_nj_ete, "Rooted Neighbor Joining Tree of HBB", "hbb_aa_nj.pdf")

# Bootstrap NJ trees
def bootstrap_nj(alignment, bs=500):
    trees = []
    n_seqs = len(alignment)
    length = alignment.get_alignment_length()
    for _ in range(bs):
        try:
            indices = np.random.choice(length, size=length, replace=True)
            sampled_aln = AlignIO.MultipleSeqAlignment(
                [SeqRecord(Seq("".join(str(rec.seq)[i] for i in indices)), id=rec.id) for rec in alignment]
            )
            dist_matrix = calculator.get_distance(sampled_aln)
            tree = constructor.nj(dist_matrix)
            trees.append(Tree(tree.format('newick'), format=1))  # Use format=1
        except Exception as e:
            print(f"Error in bootstrap iteration: {e}")
            continue
    return trees

hba_aa_nj_bs = bootstrap_nj(hba_aa_msa, bs=500)
hbb_aa_nj_bs = bootstrap_nj(hbb_aa_msa, bs=500)

# Plot bootstrapped NJ trees
def plot_bootstrap_tree(main_tree, bs_trees, title, filename):
    try:
        main_ete = Tree(main_tree.format('newick'), format=1)  # Use format=1
        for node in main_ete.traverse():
            if not node.is_leaf():
                subtree = node.get_leaf_names()
                support = sum(1 for bs_tree in bs_trees if set(subtree).issubset(set(bs_tree.get_leaf_names()))) / len(bs_trees)
                if support >= 0.5:
                    node.support = support * 100
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_support = True
        ts.title.add_face(TextFace(title, fsize=10), column=0)
        main_ete.render(filename, tree_style=ts)
    except Exception as e:
        print(f"Error plotting bootstrap tree {filename}: {e}")

plot_bootstrap_tree(hba_aa_nj, hba_aa_nj_bs, "Rooted Bootstrap NJ Tree of HBA", "hba_aa_nj_bs.pdf")
plot_bootstrap_tree(hbb_aa_nj, hbb_aa_nj_bs, "Rooted Bootstrap NJ Tree of HBB", "hbb_aa_nj_bs.pdf")

# ML trees (optional, requires PhyML)
try:
    from ete3 import PhyML
    def generate_ml_tree(alignment, output_file):
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
                SeqIO.write(alignment, tmp.name, "fasta")
                tree = PhyML()
                tree.set_alignment(tmp.name)
                tree.set_model("JTT")
                tree.set_gamma(4)
                tree.set_invariant(0.2)
                tree.optimize_branch_lengths(True)
                tree.optimize_topology(True)
                tree.run()
                ete_tree = Tree(tree.get_tree(), format=1)  # Use format=1
                return midpoint_root(ete_tree)
        except Exception as e:
            print(f"Error generating ML tree for {output_file}: {e}")
            return None
    
    hba_aa_ml = generate_ml_tree(hba_aa_msa, "hba_aa_ml.tre")
    hbb_aa_ml = generate_ml_tree(hbb_aa_msa, "hbb_aa_ml.tre")
    
    if hba_aa_ml and hbb_aa_ml:
        plot_tree(hba_aa_ml, "Rooted Maximum Likelihood Tree of HBA", "hba_aa_ml.pdf")
        plot_tree(hbb_aa_ml, "Rooted Maximum Likelihood Tree of HBB", "hbb_aa_ml.pdf")
        
        hba_aa_ml_bs = bootstrap_nj(hba_aa_msa, bs=100)
        hbb_aa_ml_bs = bootstrap_nj(hbb_aa_msa, bs=100)
        
        plot_bootstrap_tree(hba_aa_ml, hba_aa_ml_bs, "Rooted Bootstrap ML Tree of HBA", "hba_aa_ml_bs.pdf")
        plot_bootstrap_tree(hbb_aa_ml, hbb_aa_ml_bs, "Rooted Bootstrap ML Tree of HBB", "hbb_aa_ml_bs.pdf")
    else:
        print("ML tree generation failed. Using NJ trees instead.")
        hba_aa_ml = hba_aa_nj_ete
        hbb_aa_ml = hbb_aa_nj_ete
except ImportError:
    print("PhyML not installed. Skipping ML tree generation. Using NJ trees instead. Install with: conda install -c bioconda phyml")
    hba_aa_ml = hba_aa_nj_ete
    hbb_aa_ml = hbb_aa_nj_ete

# Convert scientific names to common names
def sci_to_common(names):
    common_names = []
    for name in names:
        try:
            handle = Entrez.esearch(db="taxonomy", term=name)
            record = Entrez.read(handle)
            if record["IdList"]:
                taxid = record["IdList"][0]
                handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                tax_data = Entrez.read(handle)
                common = tax_data[0].get("CommonName", name)
                common_names.append(common.replace(" ", "_"))
            else:
                common_names.append(name)
            time.sleep(0.5)  # Avoid NCBI rate limits
        except Exception as e:
            print(f"Error fetching common name for {name}: {e}")
            common_names.append(name)
    return common_names

sp1 = hba_aa_nj_ete.get_leaf_names()
sp2 = hba_aa_ml.get_leaf_names()
sp1b = hbb_aa_nj_ete.get_leaf_names()
sp2b = hbb_aa_ml.get_leaf_names()

a = sci_to_common(sp1)
a2 = sci_to_common(sp2)
b = sci_to_common(sp1b)
b2 = sci_to_common(sp2b)

print("HBA NJ common names:", a)
print("HBA ML common names:", a2)
print("HBB NJ common names:", b)
print("HBB ML common names:", b2)

# Update tree tip labels
def update_tip_labels(tree, new_labels):
    try:
        for i, leaf in enumerate(tree.get_leaves()):
            leaf.name = new_labels[i]
        return tree
    except Exception as e:
        print(f"Error updating tip labels: {e}")
        return tree

hba_aa_nj_ete = update_tip_labels(hba_aa_nj_ete, a)
hba_aa_ml = update_tip_labels(hba_aa_ml, a2)
hbb_aa_nj_ete = update_tip_labels(hbb_aa_nj_ete, b)
hbb_aa_ml = update_tip_labels(hbb_aa_ml, b2)

# Cophylogenetic plots
def plot_cophylo(tree1, tree2, title, filename):
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.rotation = 0
        tree1.render("temp1.png", tree_style=ts)
        ts.rotation = 180
        tree2.render("temp2.png", tree_style=ts)
        img1 = plt.imread("temp1.png")
        img2 = plt.imread("temp2.png")
        ax1.imshow(img1)
        ax2.imshow(img2)
        leaves1 = tree1.get_leaf_names()
        leaves2 = tree2.get_leaf_names()
        for i, label1 in enumerate(leaves1):
            for j, label2 in enumerate(leaves2):
                if label1 == label2:
                    ax1.plot([img1.shape[1], img1.shape[1] + 50], [i * img1.shape[0] / len(leaves1), i * img1.shape[0] / len(leaves1)], "k-", alpha=0.3)
                    ax2.plot([0, -50], [j * img2.shape[0] / len(leaves2), j * img2.shape[0] / len(leaves2)], "k-", alpha=0.3)
                    plt.plot([ax1.get_xlim()[1], ax2.get_xlim()[0]], [i * img1.shape[0] / len(leaves1), j * img2.shape[0] / len(leaves2)], "k-", alpha=0.3)
        ax1.set_title(f"{title} (Left: First Tree)")
        ax2.set_title(f"{title} (Right: Second Tree)")
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        os.remove("temp1.png")
        os.remove("temp2.png")
    except Exception as e:
        print(f"Error generating cophylogenetic plot {filename}: {e}")

plot_cophylo(hba_aa_nj_ete, hba_aa_ml, "NJba vs MLba", "cophylo_njba_mlba.pdf")
plot_cophylo(hbb_aa_nj_ete, hbb_aa_ml, "NJbb vs MLbb", "cophylo_njbb_mlbb.pdf")
plot_cophylo(hba_aa_ml, hbb_aa_ml, "MLba vs MLbb", "cophylo_mlba_mlbb.pdf")

# Answer questions based on tree comparisons
def compare_clades(tree1, tree2):
    try:
        clades1 = [set(node.get_leaf_names()) for node in tree1.traverse() if not node.is_leaf()]
        clades2 = [set(node.get_leaf_names()) for node in tree2.traverse() if not node.is_leaf()]
        differences = [(c1, c2) for c1 in clades1 for c2 in clades2 if c1 != c2 and c1.intersection(c2)]
        return differences
    except Exception as e:
        print(f"Error comparing clades: {e}")
        return []

# Question 1: Clade changes between HBA and HBB using ML
mlba_mlbb_diff = compare_clades(hba_aa_ml, hbb_aa_ml)
print("Clade differences between MLba and MLbb:")
for diff in mlba_mlbb_diff:
    print(f"HBA: {diff[0]}, HBB: {diff[1]}")

# Question 2: Clade changes between NJ and ML
njba_mlba_diff = compare_clades(hba_aa_nj_ete, hba_aa_ml)
njbb_mlbb_diff = compare_clades(hbb_aa_nj_ete, hbb_aa_ml)
print("Clade differences between NJba and MLba:")
for diff in njba_mlba_diff:
    print(f"NJ: {diff[0]}, ML: {diff[1]}")
print("Clade differences between NJbb and MLbb:")
for diff in njbb_mlbb_diff:
    print(f"NJ: {diff[0]}, ML: {diff[1]}")

# Question 3: Snakes in HBA trees
def check_taxonomy(tree, species, tax_group="Reptilia"):
    leaves = tree.get_leaf_names()
    reptiles = []
    for leaf in leaves:
        try:
            handle = Entrez.esearch(db="taxonomy", term=leaf.replace("_", " "))
            record = Entrez.read(handle)
            if record["IdList"]:
                taxid = record["IdList"][0]
                handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                tax_data = Entrez.read(handle)
                lineage = tax_data[0].get("Lineage", "")
                if tax_group in lineage:
                    reptiles.append(leaf)
            time.sleep(0.5)  # Avoid NCBI rate limits
        except Exception as e:
            print(f"Error checking taxonomy for {leaf}: {e}")
            continue
    is_monophyletic = any(set(reptiles).issubset(set(node.get_leaf_names())) for node in tree.traverse())
    return reptiles, is_monophyletic

hba_snakes, hba_mono = check_taxonomy(hba_aa_nj_ete, "snake")
print(f"Snakes in HBA NJ tree: {hba_snakes}, Monophyletic: {hba_mono}")
