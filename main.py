import csv
from Bio import SeqIO


# ==========================
# FUNCTIONS
# ==========================

# GC Content Function
def calculate_gc_content(sequence):

    gc_count = sequence.count("G") + sequence.count("C")

    return (gc_count / len(sequence)) * 100


# Validation Function
def validate_sequence(sequence):

    valid_nucleotides = {"A", "T", "G", "C"}

    for nucleotide in sequence:

        if nucleotide not in valid_nucleotides:

            return False

    return True


# Mutation Detection Function
def detect_mutations(seq1, seq2):

    mutations = []

    if len(seq1) != len(seq2):

        return ["Sequences have different lengths"]

    for i in range(len(seq1)):

        if seq1[i] != seq2[i]:

            mutations.append(
                f"Position {i+1}: {seq1[i]} -> {seq2[i]}"
            )

    return mutations


# ==========================
# MAIN PROGRAM
# ==========================

sequences = []

with open("results.csv", "w", newline="") as file:

    writer = csv.writer(file)

    # CSV Header
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

    # Read FASTA File
    for record in SeqIO.parse("sample.fasta", "fasta"):

        sequence = record.seq

        # Store sequence for mutation analysis
        sequences.append((record.id, str(sequence)))

        # Validation
        if validate_sequence(sequence):

            print("-" * 50)

            print("Sequence ID:", record.id)
            print("Sequence:", sequence)

            print("Length:", len(sequence))

            # GC Content
            gc_content = calculate_gc_content(sequence)

            print("GC Content:", round(gc_content, 2), "%")

            # Nucleotide Counts
            a_count = sequence.count("A")
            t_count = sequence.count("T")
            g_count = sequence.count("G")
            c_count = sequence.count("C")

            print("\nNucleotide Counts:")
            print("A:", a_count)
            print("T:", t_count)
            print("G:", g_count)
            print("C:", c_count)

            # RNA
            rna = sequence.transcribe()

            print("\nRNA Sequence:")
            print(rna)

            # Reverse Complement
            reverse_complement = sequence.reverse_complement()

            print("\nReverse Complement:")
            print(reverse_complement)

            # Protein Translation
            protein = sequence.translate()

            print("\nProtein Sequence:")
            print(protein)

            # Write Results to CSV
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

        else:

            print("-" * 50)
            print("Invalid DNA sequence detected!")
            print("Sequence ID:", record.id)


# ==========================
# MUTATION ANALYSIS
# ==========================

print("\nMutation Analysis")
print("-" * 50)

if len(sequences) >= 2:

    seq1_id, seq1 = sequences[0]
    seq2_id, seq2 = sequences[1]

    mutations = detect_mutations(seq1, seq2)

    print(f"Comparing {seq1_id} vs {seq2_id}")

    if len(mutations) > 0:

        for mutation in mutations:

            print(mutation)

    else:

        print("No mutations detected")

else:

    print("Need at least two sequences for mutation analysis")