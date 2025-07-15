from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import re

# ---------- Step 1: Load Data ----------
data_path = "/Users/921623492/Ecoli_Project/Experiment/results/combined_pangenome_matrix.csv"
df = pd.read_csv(data_path)

# ---------- Step 2: Choose gene ID ----------
gene_id = "FAHFDEJI_02221"

# ---------- Step 3: Extract and clean sequence ----------
char_Array = df[gene_id].dropna().iloc[0]
cleaned_array = [str(base) for base in char_Array if str(base).lower() in {"a", "t", "g", "c"}]
sequence = "".join(cleaned_array)

# ---------- Step 4: Format as FASTA ----------
def wrap_fasta(seq, width=70):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

fasta_query = f">{gene_id}\n{wrap_fasta(sequence)}\n"

# ---------- Step 5: Submit BLAST Query ----------
print(f"🎯 Starting BLAST for gene: {gene_id}")
result_handle = NCBIWWW.qblast(
    program="blastx",              # use 'blastn' for nucleotide, 'blastx' for protein
    database="nr",                 # 'nt' or 'refseq_protein' also work
    sequence=fasta_query,
    format_type="XML",
    hitlist_size=1
)

# ---------- Step 6: Parse Result ----------
blast_record = next(NCBIXML.parse(result_handle))

# ---------- Step 7: Function to extract gene name ----------
def extract_gene_name(description):
    # Common patterns: "...[protein/gene], gene=XYZ", or just "XYZ gene"
    match = re.search(r'\b([a-zA-Z0-9_-]+)\s+gene\b', description)
    if match:
        return match.group(1)
    
    # Try fallback patterns
    match2 = re.search(r'\bgene[:=]\s*([a-zA-Z0-9_-]+)', description)
    if match2:
        return match2.group(1)

    return None  # No gene name found

# ---------- Step 8: Print Top 5 Hits Summary ----------
print("\n✅ Top 5 BLAST Hits Summary:")
print(f"{'Rank':<5} {'Gene Name':<10} {'Identity':<12} {'e-value':<12} Description")
print("-" * 80)
for i, alignment in enumerate(blast_record.alignments[:5]):
    hsp = alignment.hsps[0]

    # Only take the first part of the description before the first '>'
    full_description = alignment.hit_def.split(">")[0].strip()

    gene_name = extract_gene_name(full_description)
    identity_pct = round(hsp.identities / hsp.align_length * 100, 1)

    print(f"{i+1:<5} {gene_name or 'None':<10} {identity_pct:<12.1f} {hsp.expect:<12.1e} {full_description}")

# ---------- Step 9: Clean Up ----------
result_handle.close()
print("\n BLAST complete.")
