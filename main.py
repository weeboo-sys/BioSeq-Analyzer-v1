import sys
import csv
from Bio import SeqIO

# Fix Windows emoji/encoding issue
sys.stdout.reconfigure(encoding='utf-8')


# ==========================
# GLOBAL MUTATION STATS
# ==========================

mutation_stats = {
    "transition": 0,
    "transversion": 0,
    "total": 0
}


# ==========================
# FUNCTIONS
# ==========================

def classify_mutation(base1, base2):
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    if base1 == base2:
        return "none"

    if (base1 in purines and base2 in purines) or (base1 in pyrimidines and base2 in pyrimidines):
        return "transition"
    else:
        return "transversion"


def detect_mutations(seq1, seq2):
    mutations = []

    if len(seq1) != len(seq2):
        return ["Sequences have different lengths"]

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:

            mutation_type = classify_mutation(seq1[i], seq2[i])

            if mutation_type == "transition":
                mutation_stats["transition"] += 1
            elif mutation_type == "transversion":
                mutation_stats["transversion"] += 1

            mutation_stats["total"] += 1

            mutations.append(
                f"Position {i+1}: {seq1[i]} -> {seq2[i]} ({mutation_type})"
            )

    return mutations


def calculate_gc_content(sequence):
    if len(sequence) == 0:
        return 0
    return (sequence.count("G") + sequence.count("C")) / len(sequence) * 100


def validate_sequence(sequence):
    valid = {"A", "T", "G", "C"}
    return all(base in valid for base in sequence)


def count_amino_acids(protein):
    counts = {}
    for aa in protein:
        counts[aa] = counts.get(aa, 0) + 1
    return counts


# ==========================
# MAIN PROGRAM
# ==========================

sequences = []

with open("results.csv", "w", newline="") as file:
    writer = csv.writer(file)

    writer.writerow([
        "Sequence ID",
        "Length",
        "GC Content (%)",
        "A Count",
        "T Count",
        "G Count",
        "C Count",
        "RNA Sequence",
        "Protein Sequence"
    ])

    for record in SeqIO.parse("sample.fasta", "fasta"):

        sequence = str(record.seq).upper()
        sequences.append((record.id, sequence))

        print("\n" + "-" * 60)

        if not validate_sequence(sequence):
            print("❌ Invalid DNA sequence detected!")
            print("Sequence ID:", record.id)
            continue

        print("🧬 Sequence ID:", record.id)
        print("Sequence:", sequence)

        print("Length:", len(sequence))

        gc_content = calculate_gc_content(sequence)
        print("GC Content:", round(gc_content, 2), "%")

        a_count = sequence.count("A")
        t_count = sequence.count("T")
        g_count = sequence.count("G")
        c_count = sequence.count("C")

        print("\nNucleotide Counts:")
        print("A:", a_count)
        print("T:", t_count)
        print("G:", g_count)
        print("C:", c_count)

        rna = record.seq.transcribe()
        protein = record.seq.translate()
        reverse_complement = record.seq.reverse_complement()

        print("\nRNA Sequence:", rna)
        print("Reverse Complement:", reverse_complement)
        print("Protein Sequence:", protein)

        amino_acid_counts = count_amino_acids(str(protein))
        print("\nAmino Acid Counts:")
        for aa, count in amino_acid_counts.items():
            print(f"{aa}: {count}")

        writer.writerow([
            record.id,
            len(sequence),
            round(gc_content, 2),
            a_count,
            t_count,
            g_count,
            c_count,
            str(rna),
            str(protein)
        ])


# ==========================
# MUTATION ANALYSIS
# ==========================

print("\n" + "=" * 60)
print("🧬 MUTATION ANALYSIS")
print("=" * 60)

if len(sequences) >= 2:

    seq1_id, seq1 = sequences[0]
    seq2_id, seq2 = sequences[1]

    print(f"Comparing {seq1_id} vs {seq2_id}\n")

    mutations = detect_mutations(seq1, seq2)

    if mutations:
        for m in mutations:
            print(m)
    else:
        print("No mutations detected")

else:
    print("Need at least two sequences for mutation analysis")


# ==========================
# MUTATION DASHBOARD
# ==========================

print("\n" + "=" * 60)
print("📊 MUTATION DASHBOARD")
print("=" * 60)

total = mutation_stats["total"]
transition = mutation_stats["transition"]
transversion = mutation_stats["transversion"]

print("Total Mutations:", total)
print("Transitions 🟢:", transition)
print("Transversions 🔴:", transversion)

if total > 0:
    print("Transition %:", round((transition / total) * 100, 2), "%")
    print("Transversion %:", round((transversion / total) * 100, 2), "%")