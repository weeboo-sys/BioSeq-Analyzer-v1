from Bio import SeqIO

for record in SeqIO.parse("sample.fasta", "fasta"):

    sequence = record.seq

    print("-" * 50)

    print("Sequence ID:", record.id)
    print("Sequence:", sequence)

    print("Length:", len(sequence))

    # GC Content
    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = (gc_count / len(sequence)) * 100

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

    print("-" * 50)