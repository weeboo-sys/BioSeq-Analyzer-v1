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


# Main Program
for record in SeqIO.parse("sample.fasta", "fasta"):

    sequence = record.seq

    # Analysis Step
    if validate_sequence(sequence):

        print("-" * 50)

        print("Sequence ID:", record.id)
        print("Sequence:", sequence)

        print("Length:", len(sequence))

        # GC Content
        gc_content = calculate_gc_content(sequence)

        print("GC Content:", round(gc_content, 2), "%")

        # Nucleotide Counts
        print("\nNucleotide Counts:")
        print("A:", sequence.count("A"))
        print("T:", sequence.count("T"))
        print("G:", sequence.count("G"))
        print("C:", sequence.count("C"))

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

    else:

        print("-" * 50)
        print("Invalid DNA sequence detected!")