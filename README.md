# BioSeq Analyzer

BioSeq Analyzer is a bioinformatics project built in Python using Biopython. The project started as a way to learn how biological sequence data can be analyzed computationally and has gradually expanded with additional sequence analysis features.

## Features

* FASTA file parsing
* DNA sequence validation
* Sequence length calculation
* GC content analysis
* Nucleotide counting (A, T, G, C)
* DNA to RNA transcription
* Protein translation
* Amino acid composition analysis
* Mutation detection between sequences
* Transition and transversion classification
* Sequence similarity calculation
* CSV export of analysis results
* Mutation statistics dashboard
* Mutation distribution visualization
* Mutation position visualization

## Technologies Used

* Python
* Biopython
* Matplotlib
* CSV
* Git
* GitHub

## Project Structure

```text
BioSeq-Analyzer/
│
├── main.py
├── sample.fasta
├── results.csv
├── mutation_distribution.png
├── mutation_positions.png
└── README.md
```

## Example Input

```fasta
>Sequence_1
ATGCGTAGCTAG

>Sequence_2
ATGAGTAGCTAG
```

## Example Output

```text
Mutation Analysis
============================================================

Sequence_1 vs Sequence_2

Position 4: C -> A (transversion)

Total Mutations: 1
Transitions: 0
Transversions: 1

Similarity: 91.67%
```

## CSV Output

The program generates a CSV file containing:

* Sequence ID
* Length
* GC Content (%)
* Nucleotide Counts
* RNA Sequence
* Protein Sequence

## What I Learned

This project helped me understand:

* Working with biological sequence data
* FASTA file handling
* Python functions and modular code
* Dictionaries and data structures
* Data export using CSV
* Mutation analysis concepts
* Data visualization with Matplotlib
* Git and GitHub workflow

## Future Improvements

* Pairwise comparison of all sequences in a FASTA file
* Similarity matrix generation
* Amino acid frequency charts
* Protein molecular weight calculation
* Integration with biological databases
* BLAST integration
* Pandas-based analysis workflow
* AI-assisted biological data analysis

## Author

Sahil Al Farish

Built as part of my journey into bioinformatics, computational biology, and AI for biological research.
