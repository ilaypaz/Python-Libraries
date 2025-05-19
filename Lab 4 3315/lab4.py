import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform, braycurtis
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from skbio.stats.ordination import pcoa
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import os
import glob
from scipy.cluster.hierarchy import leaves_list
import random
import warnings

os.chdir("/Users/ilaypaz/Downloads/Python")

# Suppress skbio RuntimeWarning about negative eigenvalues
warnings.filterwarnings("ignore", category=RuntimeWarning, module="skbio.stats.ordination")

try:
    dune = pd.read_csv("dune.csv", index_col=0)
    if dune.shape != (20, 30) or not dune.select_dtypes(include=[np.number]).equals(dune):
        raise ValueError("dune.csv must be a 20x30 matrix of numeric data")
except FileNotFoundError:
    print("Error: dune.csv not found in /Users/ilaypaz/Downloads/Python")
    print("Export from R: library(vegan); data(dune); write.csv(as.matrix(dune), 'dune.csv', row.names=TRUE)")
    exit(1)
except Exception as e:
    print(f"Error reading dune.csv: {e}")
    exit(1)

# 1. Convert distance matrices to hierarchical clustering (hclust equivalent)
dune_dist_euc = pdist(dune, metric="euclidean")
hc1 = linkage(dune_dist_euc, method="average")
dune_vd_bry = squareform([braycurtis(dune.iloc[i], dune.iloc[j]) 
                          for i in range(len(dune)) for j in range(i+1, len(dune))])
hc2 = linkage(dune_vd_bry, method="average")
dune_binary = (dune > 0).astype(int)
dune_vd_bry_bin = squareform([braycurtis(dune_binary.iloc[i], dune_binary.iloc[j]) 
                              for i in range(len(dune_binary)) for j in range(i+1, len(dune_binary))])
hc3 = linkage(dune_vd_bry_bin, method="average")

print("Euclidean hclust:\n", hc1)
print("Bray-Curtis hclust:\n", hc2)
print("Binary Bray-Curtis hclust:\n", hc3)

# 2. Convert distance matrices to matrices and save as CSV
dune_euc_matrix = squareform(dune_dist_euc)
dune_euc_df = pd.DataFrame(dune_euc_matrix, index=dune.index, columns=dune.index)
dune_euc_df.to_csv("euclidean_dist_matrix.csv")

dune_bry_df = pd.DataFrame(dune_vd_bry, index=dune.index, columns=dune.index)
dune_bry_df.to_csv("nonbinary_bray_dist_matrix.csv")

dune_bry_bin_df = pd.DataFrame(dune_vd_bry_bin, index=dune.index, columns=dune.index)
dune_bry_bin_df.to_csv("binary_bray_dist_matrix.csv")

np.random.seed(1)
dune_euc_nmds = pcoa(dune_euc_matrix)
dune_bry_nmds = pcoa(dune_vd_bry)
dune_bry_bin_nmds = pcoa(dune_vd_bry_bin)

dune_euc_points = pd.DataFrame(dune_euc_nmds.samples, index=dune.index)
dune_euc_points["Samp"] = dune_euc_points.index
dune_euc_points = dune_euc_points[["Samp"] + [col for col in dune_euc_points.columns if col != "Samp"]]

dune_bry_points = pd.DataFrame(dune_bry_nmds.samples, index=dune.index)
dune_bry_points["Samp"] = dune_bry_points.index
dune_bry_points = dune_bry_points[["Samp"] + [col for col in dune_bry_points.columns if col != "Samp"]]

dune_bry_bin_points = pd.DataFrame(dune_bry_bin_nmds.samples, index=dune.index)
dune_bry_bin_points["Samp"] = dune_bry_bin_points.index
dune_bry_bin_points = dune_bry_bin_points[["Samp"] + [col for col in dune_bry_points.columns if col != "Samp"]]

dune_euc_points.to_csv("euclidean.csv", index=False)
dune_bry_points.to_csv("nonbinaryBray.csv", index=False)
dune_bry_bin_points.to_csv("binaryBray.csv", index=False)

# 3. Convert hclust to dendrogram-like data (emulating dendro_data)
def extract_dendro_data(hc, labels):
    leaf_order = leaves_list(hc)
    segments = []
    for i, (left, right, height, _) in enumerate(hc):
        segments.append({"x": left, "y": height, "xend": right, "yend": height})
    labels_df = pd.DataFrame({
        "label": [labels[i] for i in leaf_order],
        "x": range(1, len(labels) + 1),
        "y": 0
    })
    return {"segments": pd.DataFrame(segments), "labels": labels_df}

hcd1 = extract_dendro_data(hc1, dune.index)
hcd2 = extract_dendro_data(hc2, dune.index)
hcd3 = extract_dendro_data(hc3, dune.index)

print("Euclidean dendro:\n", hcd1)
print("Bray-Curtis dendro:\n", hcd2)
print("Binary Bray-Curtis dendro:\n", hcd3)

# 4. Plot hclust and dendrogram using matplotlib
plt.figure(figsize=(8, 6))
dendrogram(hc1, labels=dune.index)
plt.title("hclust euclidean")
plt.xlabel("")
plt.ylabel("")
plt.savefig("hclust_euclidean.jpeg")
plt.close()

plt.figure(figsize=(8, 6))
dendrogram(hc2, labels=dune.index)
plt.title("hclust bray non binary")
plt.xlabel("")
plt.ylabel("")
plt.savefig("hclust_bray_nonbinary.jpeg")
plt.close()

plt.figure(figsize=(8, 6))
dendrogram(hc3, labels=dune.index)
plt.title("hclust bray binary")
plt.xlabel("")
plt.ylabel("")
plt.savefig("hclust_bray_binary.jpeg")
plt.close()

plt.figure(figsize=(8, 6))
dendrogram(hc1, labels=dune.index)
plt.title("dendrogram euclidean")
plt.xlabel("")
plt.ylabel("")
plt.savefig("dendrogram_euclidean.jpeg")
plt.close()

plt.figure(figsize=(8, 6))
dendrogram(hc2, labels=dune.index)
plt.title("dendrogram non binary bray")
plt.xlabel("")
plt.ylabel("")
plt.savefig("dendrogram_bray_nonbinary.jpeg")
plt.close()

plt.figure(figsize=(8, 6))
dendrogram(hc3, labels=dune.index)
plt.title("dendrogram binary bray")
plt.xlabel("")
plt.ylabel("")
plt.savefig("dendrogram_bray_binary.jpeg")
plt.close()

# 5. Plot NMDS using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_euc_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\neuclidean distances")
plt.tight_layout()
plt.savefig("nmds_euclidean_temp.png")
plt.close()

plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_bry_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\nBray-Curtis distances")
plt.tight_layout()
plt.savefig("nmds_bray_nonbinary_temp.png")
plt.close()

plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_bry_bin_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\nBinary Bray-Curtis distances")
plt.tight_layout()
plt.savefig("nmds_bray_binary_temp.png")
plt.close()

# 6. Export NMDS plots as PDF
plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_euc_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\neuclidean distances")
plt.tight_layout()
plt.savefig("NMDS_Dune_dataset_euclidean_distances.pdf")
plt.close()

plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_bry_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\nBray-Curtis distances")
plt.tight_layout()
plt.savefig("NMDS_Dune_dataset_Bray-Curtis_distances.pdf")
plt.close()

plt.figure(figsize=(8, 6))
sns.scatterplot(data=dune_bry_bin_points, x="PC1", y="PC2", s=50)
plt.title("NMDS Dune dataset\nBinary Bray-Curtis distances")
plt.tight_layout()
plt.savefig("NMDS_Dune_dataset_Binary_Bray-Curtis_distances.pdf")
plt.close()

# 7. Add Altitude column to labels
altitudes = ["low", "high", "sea level"]
hcd1["labels"]["Altitude"] = random.choices(altitudes, k=len(hcd1["labels"]))
hcd2["labels"]["Altitude"] = random.choices(altitudes, k=len(hcd2["labels"]))
hcd3["labels"]["Altitude"] = random.choices(altitudes, k=len(hcd3["labels"]))

# 8. Add Groups column to segments
hcd1["segments"]["Group"] = np.where(hcd1["segments"]["x"] < 10, "1", "2")
hcd2["segments"]["Group"] = np.where(hcd2["segments"]["x"] < 7, "1", "2")
hcd2["segments"]["Group"] = np.where(hcd2["segments"]["x"] > 14, "3", hcd2["segments"]["Group"])
hcd3["segments"]["Group"] = np.where(hcd3["segments"]["x"] < 7, "1", "2")
hcd3["segments"]["Group"] = np.where(hcd3["segments"]["x"] > 14, "3", hcd3["segments"]["Group"])

# 9. Plot dendro objects with colored labels and branches
def plot_dendro(hcd, title, filename):
    plt.figure(figsize=(10, 6))
    for _, row in hcd["segments"].iterrows():
        plt.plot([row["x"], row["xend"]], [row["y"], row["yend"]], 
                 color={"1": "darkred", "2": "black", "3": "darkgreen"}[row["Group"]], 
                 alpha=0.5, linewidth=1.25)
    for _, row in hcd["labels"].iterrows():
        plt.text(row["x"], row["y"] - 0.01, row["label"], rotation=90, 
                 ha="right", va="top", fontsize=8, alpha=0.5)
    plt.title(title)
    plt.xlabel("")
    plt.ylabel("")
    plt.gca().set_xticks([])
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

plot_dendro(hcd1, "Dendrogram Euclidean", "number_9_plot1.pdf")
plot_dendro(hcd2, "Dendrogram Bray-Curtis Non-Binary", "number_9_plot2.pdf")
plot_dendro(hcd3, "Dendrogram Bray-Curtis Binary", "number_9_plot3.pdf")

# 10. Export phylogenetic tree
def linkage_to_newick(hc, labels):
    tree = to_tree(hc)
    def build_newick(node, labels):
        if node.is_leaf():
            return labels[node.id]
        left = build_newick(node.left, labels)
        right = build_newick(node.right, labels)
        return f"({left},{right}):{node.dist/2}"
    return build_newick(tree, labels) + ";"

newick_str = linkage_to_newick(hc1, dune.index)
with open("prnthc.newick", "w") as f:
    f.write(newick_str)

newick_tree = Phylo.read("prnthc.newick", "newick")
plt.figure(figsize=(8, 6))
Phylo.draw(newick_tree, do_show=False)
plt.savefig("newick_tree.jpeg")
plt.close()

with open("nxs.nexus", "w") as f:
    f.write("#NEXUS\nbegin trees;\ntree tree1 = " + newick_str + "\nend;")
nexus_tree = Phylo.read("nxs.nexus", "nexus")
plt.figure(figsize=(8, 6))
Phylo.draw(nexus_tree, do_show=False)
plt.savefig("nexus_tree.jpeg")
plt.close()

# 11. csvreader function
def csvreader(folder_path):
    files = glob.glob(os.path.join(folder_path, "*.csv"))
    matrices = []
    skip_files = ['euclidean.csv', 'nonbinaryBray.csv', 'binaryBray.csv']
    for file in files:
        if os.path.basename(file) in skip_files:
            print(f"Skipping {file}: Contains NMDS coordinates, not raw data")
            continue
        try:
            data = pd.read_csv(file)
            if 'Samp' in data.columns:
                data = data.drop(columns=['Samp'])
            numeric_data = data.select_dtypes(include=[np.number])
            if numeric_data.empty:
                print(f"Skipping {file}: No numeric data")
                continue
            log_data = np.log1p(numeric_data)
            dm = squareform([braycurtis(log_data.iloc[i], log_data.iloc[j]) 
                             for i in range(len(log_data)) for j in range(i+1, len(log_data))])
            if np.any(np.isnan(dm)) or not np.allclose(dm, dm.T):
                print(f"Skipping {file}: Data must be symmetric and cannot contain NaNs")
                continue
            nmds = pcoa(dm)
            # Calculate stress (goodness-of-fit)
            stress = np.sqrt(((nmds.samples.values - nmds.samples.values.mean(axis=0))**2).sum() / 
                             ((dm - dm.mean())**2).sum())
            print(f"Processed {file}: NMDS stress = {stress:.4f}")
            matrices.append(nmds)
        except Exception as e:
            print(f"Error processing {file}: {e}")
    return matrices

d = csvreader("/Users/ilaypaz/Downloads/Python")
print("NMDS results:", len(d), "matrices")
