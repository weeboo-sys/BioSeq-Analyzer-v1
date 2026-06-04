import csv
from Bio import SeqIO
import matplotlib.pyplot as plt


# ==========================
# MUTATION STATS
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

            mtype = classify_mutation(seq1[i], seq2[i])

            if mtype == "transition":
                mutation_stats["transition"] += 1
            elif mtype == "transversion":
                mutation_stats["transversion"] += 1

            mutation_stats["total"] += 1

            mutations.append(f"Position {i+1}: {seq1[i]} -> {seq2[i]} ({mtype})")

    return mutations


def calculate_gc_content(seq):
    if len(seq) == 0:
        return 0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def validate_sequence(seq):
    return all(base in {"A", "T", "G", "C"} for base in seq)


def count_amino_acids(protein):
    counts = {}
    for aa in protein:
        counts[aa] = counts.get(aa, 0) + 1
    return counts


# ==========================
# GRAPH FUNCTIONS
# ==========================

def plot_mutation_stats(stats):
    labels = ["Transitions", "Transversions"]
    values = [stats["transition"], stats["transversion"]]

    plt.figure()
    plt.bar(labels, values)
    plt.title("Mutation Type Distribution")
    plt.ylabel("Count")
    plt.show()


def plot_mutation_positions(seq1, seq2):
    positions = []

    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            positions.append(i)

    plt.figure()
    plt.scatter(positions, [1] * len(positions))
    plt.title("Mutation Positions Along Sequence")
    plt.yticks([])
    plt.xlabel("Position in Sequence")
    plt.show()


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

        seq = str(record.seq).upper()
        sequences.append((record.id, seq))

        print("\n" + "-" * 60)
        print(" ID:", record.id)

        if not validate_sequence(seq):
            print(" Invalid sequence")
            continue

        gc = calculate_gc_content(seq)
        print("GC%:", round(gc, 2))

        rna = record.seq.transcribe()
        protein = record.seq.translate()

        print("RNA:", rna)
        print("Protein:", protein)

        writer.writerow([
            record.id,
            len(seq),
            round(gc, 2),
            seq.count("A"),
            seq.count("T"),
            seq.count("G"),
            seq.count("C"),
            str(rna),
            str(protein)
        ])


# ==========================
# MUTATION ANALYSIS
# ==========================

print("\n" + "=" * 60)
print(" MUTATION ANALYSIS")
print("=" * 60)

if len(sequences) >= 2:

    s1_id, s1 = sequences[0]
    s2_id, s2 = sequences[1]

    print(f"{s1_id} vs {s2_id}")

    mutations = detect_mutations(s1, s2)

    for m in mutations:
        print(m)

    print("\n DASHBOARD")
    total = mutation_stats["total"]

    print("Total:", total)
    print("Transitions :", mutation_stats["transition"])
    print("Transversions :", mutation_stats["transversion"])

    if total > 0:
        print("Transition %:", round(mutation_stats["transition"] / total * 100, 2))
        print("Transversion %:", round(mutation_stats["transversion"] / total * 100, 2))

    # ==========================
    # GRAPHS
    # ==========================

    plot_mutation_stats(mutation_stats)
    plot_mutation_positions(s1, s2)

else:
    print("Need at least two sequences for mutation analysis")