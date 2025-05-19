from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import urllib.request
import random
import os
from collections import defaultdict

# 1. Function to read FASTA files and count sequences and their lengths
def read_fasta_count_sequences(fasta_file):
    sequences = []
    # Download the FASTA file if it's a URL
    if fasta_file.startswith("http"):
        local_file = "practice_nuc.fa"
        urllib.request.urlretrieve(fasta_file, local_file)
        fasta_file = local_file
    
    # Read FASTA file and count sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))
    
    # Count sequences and lengths
    seq_count = len(sequences)
    seq_lengths = [len(seq[1]) for seq in sequences]
    
    return seq_count, sequences, seq_lengths

# Example usage with the provided FASTA file
fasta_url = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/practice_nuc.fa"
seq_count, sequences, seq_lengths = read_fasta_count_sequences(fasta_url)
print(f"Number of sequences: {seq_count}")
print(f"Sequence lengths: {seq_lengths}")
print(f"Sequences: {sequences}")

# 2. Function to reverse complement a nucleotide sequence (DNA or RNA)
def reverse_complement(seq, is_rna=False):
    # Convert to Seq object
    seq = Seq(seq)
    # Get reverse complement
    rev_comp = seq.reverse_complement()
    # Convert to string for further processing
    rev_comp_str = str(rev_comp)
    # Ensure output is DNA (replace U with T if RNA)
    if is_rna:
        rev_comp_str = rev_comp_str.replace("U", "T")
    return rev_comp_str

# Example sequences from the FASTA file
Dseq = "".join(sequences[0][1])  # First sequence as DNA
Rseq = Dseq.replace("T", "U")     # Convert to RNA
print(f"Reverse complement (DNA): {reverse_complement(Dseq, is_rna=False)}")
print(f"Reverse complement (RNA): {reverse_complement(Rseq, is_rna=True)}")

# 3. Generate five 500-amino-acid strings using the standard 20 AA alphabet
AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")  # Standard 20 amino acids
my_fave_proteins = [""] * 5
for i in range(5):
    # Ensure non-identical sequences by sampling randomly
    random.seed(i)  # Set seed for reproducibility and variation
    my_fave_proteins[i] = "".join(random.choices(AA_ALPHABET, k=500))
print(my_fave_proteins)

# 4. Assign my_fave_proteins to a SeqRecord list (equivalent to XStringSet)
my_fave_proteins_set = [
    SeqRecord(Seq(seq), id=f"protein_{i+1}", annotations={"molecule_type": "protein"})
    for i, seq in enumerate(my_fave_proteins)
]
print(my_fave_proteins_set)

# 5. Change the first AA to 'M' and assign creative names
creative_names = ["Beelzebub", "Lucifer", "Belial", "Behemoth", "Asmodeus"]
for i, record in enumerate(my_fave_proteins_set):
    record.id = creative_names[i]
    record.name = creative_names[i]
    record.description = creative_names[i]
    # Replace first AA with 'M'
    seq = str(record.seq)
    record.seq = Seq("M" + seq[1:])
print(my_fave_proteins_set)

# 6. Get AA frequency and plot a bar plot
def aa_frequency(seq_records):
    freqs = []
    for record in seq_records:
        seq = str(record.seq)
        # Count frequency of each AA
        aa_counts = defaultdict(int)
        for aa in seq:
            aa_counts[aa] += 1
        # Convert to proportions
        total = len(seq)
        aa_freq = {aa: count / total * 100 for aa, count in aa_counts.items()}
        aa_freq["Sqs"] = record.id
        freqs.append(aa_freq)
    return pd.DataFrame(freqs)

# Calculate frequencies
freq_df = aa_frequency(my_fave_proteins_set)
# Melt dataframe to long format (like tidyr::gather)
freq_df = freq_df.melt(id_vars=["Sqs"], var_name="Bases", value_name="Frequency")
# Filter to standard AAs (if needed) and remove NaN
freq_df = freq_df.dropna()

# Plot bar plot using seaborn
plt.figure(figsize=(12, 8))
sns.barplot(x="Bases", y="Frequency", hue="Sqs", data=freq_df)
plt.title("Amino Acid Frequencies Across Proteins")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency (%)")
plt.legend(title="Proteins")
plt.tight_layout()
plot = plt.gcf()
plt.show()

# Save plot as PDF
plot.savefig("aa_frequencies.pdf", format="pdf", bbox_inches="tight", dpi=300)
plt.close()

# 7. Export my_fave_proteins_set as a FASTA file
SeqIO.write(my_fave_proteins_set, "ilay_fave_prots.fasta", "fasta")

# 8. Import first 200 Archaea sequences from GenBank
Entrez.email = "pazlilay@gmail.com"  # Replace with your actual email
Entrez.api_key = "3102925945433d6a3dd2b44772a890ea8a08"          # Replace with your NCBI API key
def fetch_archaea_sequences(max_seqs=50):  # Reduced to 50 to avoid rate limits
    try:
        search_query = "Archaea[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_seqs)
        record = Entrez.read(handle)
        handle.close()
        seq_ids = record["IdList"]
        total_seqs = int(record["Count"])
        
        # Fetch sequences
        sequences = []
        handle = Entrez.efetch(db="nucleotide", id=seq_ids, rettype="fasta", retmode="text")
        for seq_record in SeqIO.parse(handle, "fasta"):
            sequences.append(seq_record)
        handle.close()
        
        return sequences, total_seqs
    except Exception as e:
        print(f"Error fetching GenBank sequences: {e}")
        return [], 0

archaea_seqs, total_archaea = fetch_archaea_sequences()
print(f"Total Archaeal sequences in GenBank: {total_archaea}")
print(f"Fetched sequences: {len(archaea_seqs)}")

# 9. Assign to SeqRecord list and reverse complement
my_fave_arch_set = []
for seq_record in archaea_seqs:
    # Reverse complement
    rev_comp_seq = seq_record.seq.reverse_complement()
    new_record = SeqRecord(
        rev_comp_seq,
        id=seq_record.id,
        name=seq_record.name,
        description=seq_record.description,
        annotations={"molecule_type": "DNA"}
    )
    my_fave_arch_set.append(new_record)
print(my_fave_arch_set)

# 10. Find the length of the longest sequence
longest_length = max(len(record.seq) for record in my_fave_arch_set) if my_fave_arch_set else 0
print(f"Length of the longest sequence: {longest_length}")

# 11. Export my_fave_arch_set as a FASTA file
if my_fave_arch_set:
    SeqIO.write(my_fave_arch_set, "ilay_fave_arch.fa", "fasta")
else:
    print("No archaeal sequences to export")
