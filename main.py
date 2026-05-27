import csv
from Bio import SeqIO


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


# Open CSV File
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
        "C Count"
    ])

    # Main Program
    for record in SeqIO.parse("sample.fasta", "fasta"):

        sequence = record.seq

        # Validation Check
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

            # RNA Conversion
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

            # Write To CSV
            writer.writerow([
                record.id,
                len(sequence),
                round(gc_content, 2),
                a_count,
                t_count,
                g_count,
                c_count
            ])

        else:

            print("-" * 50)
            print("Invalid DNA sequence detected!")