from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

filepath = "viruses/data"

# Nuskaitykite FASTA failus
def read_fasta(filepath):
    sequences = []
    for record in SeqIO.parse(filepath, "fasta"):
        sequences.append(record.seq)
    return sequences

# Funkcija surasti start ir stop kodonus
def find_start_stop_pairs(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    found_pairs = []
    
    i = 0
    while i < len(sequence) - 2:
        codon = sequence[i:i+3]
        if codon == start_codon:
            for j in range(i + 3, len(sequence) - 2, 3):
                next_codon = sequence[j:j+3]
                if next_codon in stop_codons:
                    found_pairs.append((i, j + 3))
                    i = j
                    break
        i += 3
    return found_pairs

# Reverse komplemento sekos suradimas
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

# Filtruoti fragmentus trumpesnius nei 100 bp
def filter_short_sequences(pairs, min_length=100):
    return [(start, stop) for start, stop in pairs if stop - start >= min_length]

# DNR sekos konvertavimas į baltymus
def translate_to_protein(sequence, pairs):
    proteins = []
    for start, stop in pairs:
        coding_sequence = sequence[start:stop]
        protein = Seq(coding_sequence).translate()
        proteins.append(str(protein))
    return proteins

# Kodonų ir dikodonų dažnio analizė
def codon_dicodon_frequency(protein_sequences):
    codon_count = defaultdict(int)
    dicodon_count = defaultdict(int)
    
    for protein in protein_sequences:
        for i in range(0, len(protein), 1):
            codon = protein[i:i+1]
            codon_count[codon] += 1
            if i + 1 < len(protein):
                dicodon = protein[i:i+2]
                dicodon_count[dicodon] += 1
    
    return codon_count, dicodon_count

# Atstumo matricos sudarymas (pvz., naudojant Euklido atstumą)
def calculate_distance_matrix(frequencies):
    import numpy as np
    matrix = np.zeros((len(frequencies), len(frequencies)))
    for i, freq1 in enumerate(frequencies):
        for j, freq2 in enumerate(frequencies):
            if i != j:
                dist = np.linalg.norm(np.array(freq1) - np.array(freq2))
                matrix[i][j] = dist
    return matrix

# Pagrindinis funkcijų vykdymas
fasta_file = "virus_sequences.fasta"
sequences = read_fasta(fasta_file)

# Analizuoti kiekvieną seką
for sequence in sequences:
    reverse_seq = reverse_complement(sequence)
    
    pairs = find_start_stop_pairs(sequence)
    reverse_pairs = find_start_stop_pairs(reverse_seq)
    
    filtered_pairs = filter_short_sequences(pairs)
    filtered_reverse_pairs = filter_short_sequences(reverse_pairs)
    
    protein_sequences = translate_to_protein(sequence, filtered_pairs)
    reverse_protein_sequences = translate_to_protein(reverse_seq, filtered_reverse_pairs)
    
    codon_freq, dicodon_freq = codon_dicodon_frequency(protein_sequences + reverse_protein_sequences)
    
    # Pavyzdžiui, čia išveskite kodonų dažnius
    print("Kodonų dažnis:", codon_freq)
    print("Dikodonų dažnis:", dicodon_freq)

# Sudaryti atstumo matricą (galite pritaikyti klasterizavimo metodą)
dist_matrix = calculate_distance_matrix([codon_freq, dicodon_freq])
