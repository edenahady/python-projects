'''Eden Ahady 318948106'''
from sys import argv


class Cell:
    def __init__(self, name, genome, reading_frames):
        self.name = name
        self.genome = genome
        self.reading_frames = reading_frames
    #finding the repeat sequences in the genome
    def get_srr(self,i):
        """
             :param self: our cell
             :param i : the index of the genome
             :return: string representing
             """
        index = i % len(self.genome) #makes the index circular as asked
        dna = self.genome[index]
        dict = {}
        # loop over all the possible sizes of the srr
        for length in range(1, 7):
            # check for all index
            for i in range(len(dna)):
                temp_i = i + length
                srr = dna[i:temp_i]
                times = 1
                # if the sequence is equal anf if it's repeating
                while dna[temp_i:temp_i + length] == srr and temp_i < len(dna):
                    temp_i = temp_i + length
                    times = times + 1
                    # put the sequence in the dictionary if apear more than 3 times
                    if times >= 3 and srr not in dict:
                        dict[srr] = times
                    elif srr in dict:
                        if times > dict[srr]:
                            dict[srr] = times
        # sort the dictionary by the number of repeats
        srr_list = sorted(dict.items(), key=lambda x: x[1])
        return srr_list
    #dna into rna
    def transcribe(self, i):
        """
                     :param self: our cell
                     :param i : the index of the genome
                     :return: rna
                     """
        dna = self.genome[i]
        rna = ""
        upper_dna = dna.upper()  # makes sure seq is all uppercase
        # replace all bases
        for i in range(len(upper_dna)):
            if upper_dna[i] == 'A':
                rna += 'U'
            elif upper_dna[i] == 'T':
                rna += 'A'
            elif upper_dna[i] == 'C':
                rna += 'G'
            elif upper_dna[i] == 'G':
                rna += 'C'
        # return the rna sequence
        return rna[::-1]

    def get_protein(seq):
        """
              :param seq: our sequence
              :return: the translated protein
              """
        # table of proteins according to their codons
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
        # start codon
        codon_start = "AUG"
        # list of 3 stop codons
        list_stop = ["UAG", "UGA", "UAA"]
        # find the first AUG
        start_op = seq.find(codon_start)
        longest_seq = ""
        # there is no AUG
        if start_op == -1:
            return None
        while start_op >= 0:
            temp_seq = "AUG"
            # starts from the first AUG till the end
            seq = seq[start_op + 3:]
            i = 0
            while (i < len(seq)):
                # the next codon
                codon = seq[i:i + 3]
                # if we reached stop codon or to the end
                if codon in list_stop or len(codon) < 3:
                    # replacing with the longest sequence
                    if len(longest_seq) < len(temp_seq):
                        longest_seq = temp_seq
                        break
                else:
                    temp_seq += codon
                i += 3
            # replacing with the longest sequence
            if len(longest_seq) < len(temp_seq):
                longest_seq = temp_seq
            else:
                break
            # finding the next AUG again
            start_op = seq.find(codon_start)
        protein = ""
        # replace the codon with the right letter
        for i in range(0, len(longest_seq), 3):
            codon = longest_seq[i:i + 3]
            if i+3 >= len(longest_seq):
                protein += table[codon]
            else:
                protein += table[codon] + ";"
        return protein

    def translate(self,i):
        """
               :param seq: the rna sequence we want to translate
               :param i: the index of the genome
               :return: the biggest open frame that can be translated
               """
        frame = self.genome[i]
        seq = self.transcribe(i)
        k = int(self.reading_frames[i])
        # cutting the bases acording to the open frame given
        if (int(k) % 3 == 1):
            pro1 = Cell.get_protein(seq)
            return pro1

        if (int(k) % 3 == 2):
            seq = seq[1:len(seq) - 2]
            pro2 = Cell.get_protein(seq)
            return pro2

        if (int(k) % 3 == 0):
            seq = seq[2:len(seq) - 1]
            pro3 = Cell.get_protein(seq)
            return pro3

    def repertoire(self):
        """
            :param self: is our cell
            :return: list of tuples with the srr protein
        """
        list = []
        for i in range(len(self.genome)): #check for every sequence in the genome
            srr = self.get_srr(i)
            protien = self.translate(i)
            list.append((srr, protien))
        return list

class StemCell(Cell):
    def __init__(self, my_genome, reading_frames):
        super().__init__(name="Stem Cell" , genome=my_genome, reading_frames=reading_frames)
    #makes s deep copy of the cell
    def copy(self):
        """
         :param self: is our stem cell
        :return: the deep copy of the cell
        """
        dup = Cell(self.name, self.genome, self.reading_frames)
        return dup

    def __mul__(self, coafficient):
        """

        :param coafficient: is the coafficient of the signal
        :return: overrides the multiplication function
        """
        cells_list = [self]
        for i in range(coafficient):
            duplicate = self.copy()
            cells_list.append(duplicate)

    def differentiate(self, name, params):
        """

        :param name: is the genome
        :param params: list of parameters from the file
        :return: stem cell into nerve cell or muscle cell
        """
        our_cell = None
        if len(params) == 1 and name == "NC":
            our_cell = self.NerveCell(float(params[0]), self) #nerve cell
        if len(params) == 2 and name == "MC":
            our_cell = self.MuscleCell(self, params[0], float(params[1]))#muscle cell

        return our_cell

    def mitosis(self, cell):
        """

        :param cell: is our cell
        :return: making mitosis by makin a list with cell the function called and a copy
        """
        dupe = Cell(cell.name, cell.genome, cell.reading_frames)
        cells_list = [cell, dupe]
        return cells_list

    class NerveCell(Cell):
        def __init__(self, coefficient: float, stemcell):
            super().__init__(name="Nerve Cell", genome=stemcell.genome, reading_frames=stemcell.reading_frames)
            self.coefficient = coefficient
            self.rec_signal = None

        def recieve(self, our_signal):
            """

            :param our_signal: is the given signal
            """
            self.rec_signal = our_signal

        def send(self, signal):
            """

            :param signal: is the given signal
            :return: return the signal*it's coafficient
            """
            return signal * self.coefficient

    class MuscleCell(Cell):
        def __init__(self, stemcell, file, threshold):
            super().__init__(name="Muscle Cell", genome=stemcell.genome, reading_frames=stemcell.reading_frames)
            self.file = file
            self.threshold = threshold
            self.signal = None

        def recieve(self, our_signal):
            """

            :param our_signal: is the given signal
            :return:if the signal is bigger than the threshold it writes to a file
            """
            if our_signal > self.threshold:
                with open(self.file, 'a') as fi:
                    fi.write(str(our_signal) + ", " + "I like to move it\n")
                self.signal = our_signal

class NerveNetwork:
    def __init__(self, nerve_list, muscle):
        self.nerve_list = nerve_list
        self.muscle = muscle

    def send_signal(self, our_signal):
        """

        :param our_signal: is the given signal
        :return: send the signal from nerve cell to another and them to muscle cell
        """
        for c in self.nerve_list:
            c.recieve(our_signal)
            our_signal = c.send(our_signal) #send the signal
        self.muscle.recieve(our_signal)

def print_srr(tup_list):
    """

    :param tup_list: is the list of tuples
    :return: printing the srrs
    """
    if len(tup_list) > 0: #list is not empty
        for i in range(len(tup_list)):
            print(','.join(str(x) for x in tup_list[i]), sep='', end='') #prints
            if i < len(tup_list)-1: #makes sure not to print ; after last protien
                print(";", sep='', end='')
        print('')
    else: #list is empty
        print("No simple repeats in DNA sequence")

def assert_check(divided_line, divided_dna, divided_rf , divided_param):
    """

    :param divided_line: is the list of lines
    :param divided_dna: is the divided dna column
    :param divided_rf: is the divided reading frames column
    :param divided_param: is the divided parameters column
    """
    #if there are 4 colums
    assert len(divided_line) == 4, "File illegal"
    #check the name of the cell
    assert divided_line[0] == 'MC' or divided_line[0] == 'NC', "File illegal"
    for dna in divided_dna:
        for letter in range(len(dna)):
            #check if all the sequence contains only ATGC
            assert dna[letter] == 'A' or dna[letter] == 'T' or dna[letter] == 'C' or dna[letter] == 'G',\
                "File illegal"
    #length og the reading frames equal to dna seq
    assert len(divided_rf) == len(divided_dna), "File illegal"
    for rf in divided_rf:
        #if the number has no decimal point
        assert rf in ['1', '2', '3'], "File illegal"
    #split_param = divided_line[3].split(',')
    #int(divided_param[1]))) or
    assert divided_line[0] == "MC" and len(divided_param) == 2 and (float(divided_param[1]) > 0) or \
           (divided_line[0] == "NC" and len(divided_param) == 1 and float(divided_param[0]) > 0), "File illegal"

def main(argv):
    """

    :param argv: is the list of arguments
    """
    nerve_list = []
    with open(argv[1], 'r') as tsb_file:
        lines = tsb_file.readlines()#open the file
        for line in lines[1:]:
            divided_line = line.split('\t') #split by tabs
            divided_dna = divided_line[1].split(',') #split the dna columns
            divided_rf = divided_line[2].split(',') #split the dna reading frames
            divided_param = divided_line[3].split(',') #split parameters columns
            assert_check(divided_line, divided_dna, divided_rf,divided_param)
            stemcell = StemCell(divided_dna, divided_rf)
            if divided_line[0] != "NC":
                muscleCell = stemcell.differentiate("MC", divided_param)#stem cell into muscle
                reps = stemcell.repertoire()#doing repertoire
                for rep in reps:
                    tup_list = rep[0]
                    print_srr(tup_list) #print
                    if rep[1] is not None:
                        print("Translation:", rep[1])
                    else:
                        print("Non-coding RNA") #printing the protein
            else:
                    cells = stemcell.mitosis(stemcell) #doing mitosis for stem cell
                    for cell in cells:
                        cell_sort = stemcell.differentiate("NC", divided_param)
                        nerve_list.append(cell_sort)

        nervenet = NerveNetwork(nerve_list, muscleCell)
        split_arg = argv[2].split(',')
        for n in split_arg:
            signal = float(n)
            nervenet.send_signal(signal)
if __name__ == '__main__':
    main(argv)