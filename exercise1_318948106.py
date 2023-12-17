def reverse_transcribe(rna):
    dna = rna
    dna_upp = dna.upper()
    dna_upp = dna_upp.replace("A", "T")
    dna_upp = dna_upp.replace("U", "A")
    dna_upp = dna_upp.replace("C", "G")
    dna_upp = dna_upp.replace("G", "C")
    return dna_upp[::-1]

def read_frame(seq, frame):
    if (frame % 3 == 1):
       pro1 = translate(seq)
       return pro1

    if(frame % 3 == 2):
        seq = seq[1:len(seq)-2]
        pro2 = translate(seq)
        return pro2

    if (frame % 3 == 0):
        seq = seq[2:len(seq) - 1]
        pro3 = translate(seq)
        return pro3


def translate(seq):
    table = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '', 'UAG': '',
        'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]

    return protein
def main(argv):
    print("DNA sequence:", reverse_transcribe(argv[1]))
    if translate(argv[2], argv[3]) is None:
        print("Non-coding RNA")
    else:
        print(translate(argv[2], argv[3]))


if _name_ == '_main_':
    main(argv)
