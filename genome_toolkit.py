from Bio import SeqIO

class GenomeToolkit:
    """A toolkit for analysing DNA/RNA sequences."""

    def __init__(self):
        "Initialise attributes for the toolkit."
        print("------ Genome Toolkit ------")
        print("\n - Developed by Aaron Ngai")
        print(" - https://github.com/angai4")
        print("\n----------------------------\n")

    def dna_to_rna(self, sequence):
        """
        This method takes a DNA sequence corresponding to a coding strand and
        returns the RNA sequence
        """
        rna_seq = sequence.upper().replace('T', 'U')
        return rna_seq
    
    def count_nts(self, sequence):
        """
        This method counts the occurences of each nucleotide in a sequence
        """
        sequence = sequence.upper()
        A, C, G, T = 0, 0, 0, 0
        for nt in sequence:
            if nt == 'A':
                A += 1
            elif nt == 'C':
                C += 1
            elif nt == 'G':
                G += 1                
            elif nt == 'T':
                T += 1
        print(f"A: {A}")
        print(f"C: {C}")  
        print(f"G: {G}")  
        print(f"T: {T}")

    def rev_comp(self, sequence):
        """
        This method takes in a DNA sequence and returns its reverse complement 
        """                  
        comp = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C'
        }
        sequence = sequence.upper()
        rev_comp_seq = ""
        for nt in reversed(sequence):
            rev_comp_seq += comp[nt]
        return rev_comp_seq
    
    def translation(self, sequence):
        """
        This method translates a nucleic acid sequence into a protein sequnce, 
        until the end or until it comes across a stop codon
        """
        genetic_code = {
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',    # Serine
            'UUC': 'F', 'UUU': 'F',    # Phenylalanine
            'UUA': 'L', 'UUG': 'L',    # Leucine
            'UAC': 'Y', 'UAU': 'Y',    # Tyrosine
            'UAA': None, 'UAG': None,    # Stop
            'UGC': 'C', 'UGU': 'C',    # Cysteine
            'UGA': None,    # Stop
            'UGG': 'W',    # Tryptophan
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',    # Leucine
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',    # Proline
            'CAC': 'H', 'CAU': 'H',    # Histidine
            'CAA': 'Q', 'CAG': 'Q',    # Glutamine
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',    # Arginine
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I',    # Isoleucine
            'AUG': 'M',    # Methionine (Start)
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',    # Threonine
            'AAC': 'N', 'AAU': 'N',    # Asparagine
            'AAA': 'K', 'AAG': 'K',    # Lysine
            'AGC': 'S', 'AGU': 'S',    # Serine
            'AGA': 'R', 'AGG': 'R',    # Arginine
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',    # Valine
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',    # Alanine
            'GAC': 'D', 'GAU': 'D',    # Aspartic Acid
            'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G'     # Glycine
        }
        sequence = sequence.upper().replace('T', 'U')
        protein_seq = ""
        i = 0
        while i+2 < len(sequence):
            codon = sequence[i:i+3]
            amino_acid = genetic_code[codon]
            if amino_acid is None:
                break
            protein_seq = protein_seq + amino_acid
            i += 3
        return protein_seq
    
    def calc_gc(self, fasta_file):
        """
        This function reads a fasta file as input and prints the sequence with 
        the highest GC% aswell as the sequence ID
        """
        file = open(fasta_file, 'r')
        records = {}
        for record in SeqIO.parse(file, 'fasta'):
            records[record.id] = str(record.seq)
        file.close()
        scores = []
        score_seq_pair = {}
        for seqs in records.values():
            num_gc = seqs.count('G') + seqs.count('C')
            value = num_gc / len(seqs) * 100
            scores.append(value)
            score_seq_pair[value] = seqs
        max_gc_seq = score_seq_pair[max(scores)]
        target_v = max_gc_seq
        for k, v in record.items():
            if v == target_v:
                k_for_v = k
                break
        print(f"Highest GC%: {round(max(scores), 3)}%")
        print(f"Sequence ID: {k_for_v}")