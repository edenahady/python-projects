#Eden Ahady 318948106
import re
from sys import argv
import subprocess
from matplotlib import pyplot as plt
from Bio import SeqIO


def translate_to_regex(pro_pattern):
    """
    :param pro_list: is the prosite patterns
    Loads the PROSITE patterns from the text file and convert them to regex
    :return: regex list that matches the prosite pattern
    """
    #replace { with [^
    regex = pro_pattern.replace("{", "[^")
    #replace } with ]
    regex = regex.replace("}", "]")
    #replace ( with {
    regex = regex.replace("(", "{")
    # replace ) with }
    regex = regex.replace(")", "}")
    # remove -
    regex = regex.replace("-", "")
    # replace x with .
    regex = regex.replace("x", ".")
    # replace < with ^
    regex = regex.replace("<", "^")
    # replace > with $
    regex = regex.replace(">", "$")
    # replace PROSITE character class notation with regular expression syntax
    return regex


def count_occurrences(regex, fasta_seq):
    """
    :param regex: is the rexes we want to count it's occurrences
    :param fasta_seq: is the fasta we look in
    :return: Counts the occurrences of each regular expression in the FASTA sequences
    """
    occurrences = []
    hits = []
    # find all matches
    for record in SeqIO.parse(fasta_seq, "fasta"):
        hits.append(str(record.seq))
    for hit in hits:
        flag = re.findall(regex, hit)
        # if we have a match
        if flag:
            # add matches to list
            occurrences.extend(flag)
    return occurrences

def occurrences_length(second_fasta, third_fasta, regex):
    """
    :param second_fasta: the second fasta file
    :param third_fasta: the third fasta file
    :param regex: is the regex we look for
    :return: the length of the occurrences and how many occurrences for each length
    """
    dict_2 = {}
    dict_3 = {}
    # check the number of occurrences for every motive
    for r in regex:
        count = count_occurrences(r, second_fasta)
        dict_2 = {**dict_2, **get_length(count)}
        count = count_occurrences(r, third_fasta)
        dict_3 = {**dict_3, **get_length(count)}
    return dict_2, dict_3

def get_length(matches):
    """
    :param matches: list of matches
    :return: dictionary with length of matches as keys and occurrences as values
    """
    lengths = {}
    for match in matches:
        length = len(match)
        # we add 1 to every cell in the dict that represents the pattern
        lengths[length] = lengths.get(length, 0) + 1
    return lengths


def plot_occurrences(dict, dict_2, dict_3):
    """
    :param dict: is the occurrence's dictionary
    :param dict_2: is the second fasta dictionary
    :param dict_3: is the third fasta dictionary
     Plots a bar graph comparing the number of occurrences of the PROSITE templates
    """
    f, axis = plt.subplots(2)
    f.set_figheight(10)
    f.set_figwidth(10)

    # upper plot
    axis[0].set_ylabel("Number of matches")
    axis[0].set_xlabel("Pattern name")
    axis[0].set_title("Number of matches in the prosite sequence")
    axis[0].bar(dict.keys(), dict.values())

    # lower plot
    axis[1].scatter(dict_2.keys(), dict_2.values(), c="blue", label="fasta file 1")
    axis[1].scatter(dict_3.keys(), dict_3.values(), c="red", label="fasta file 2")
    axis[1].set_ylabel("Quantity")
    axis[1].set_xlabel("length of pattern")
    axis[1].set_title("length of pattern in the fasta and their quantity")
    axis[1].set_yscale('log')
    axis[1].legend()
    plt.savefig("318948106.png")

def subprocessing(path_fastp, fastq):
    """

    :param path_fastp: us the fastp path
    :param fastq: is our fastq file
    runs the program on the file givven
    """
    before_change = ""
    after_change = ""
    result = subprocess.run([path_fastp, "--in1", fastq, "--out1", "/dev/null", "-j", "/dev/null", "-h", "/dev/null"],
                         capture_output=True, text=True, timeout=5)
    if(result.returncode != 0):
        raise Exception(result.stderr)
    split_lines = str(result.stderr).split("\n")
    output = [line.split(": ")[1] for line in split_lines if "total reads" in line]
    if len(output) == 2:
        before_change, after_change = output[:2]
    else:
        raise Exception(result.stderr)
    # check if the change is valif and calculate the diff
    if before_change.isdigit() and after_change.isdigit():
        return int(before_change) - int(after_change)
    else:
        raise Exception(result.stderr)


def main(argv):
    dict = {}
    regexs = []
    patterns = []
    #reading the lines of the file
    with open(argv[1]) as f:
        list_lines = f.read().splitlines()  # list of lines
        first_fasta = list_lines[0]
        second_fasta = list_lines[1]
        third_fasta = list_lines[2]
        fastq = list_lines[3]
        #split the lines by the delimeter
        for line in list_lines[4:]:
            split_dlimeter = tuple(line.split("; "))
            patterns.append(split_dlimeter)
    for pattern, p in patterns:
        regex = translate_to_regex(pattern)
        regexs.append(regex)
        count = count_occurrences(regex, first_fasta)
        dict[p] = len(count)
    # we count motives in the other fastas files
    dict_2, dict_3 = occurrences_length(second_fasta, third_fasta, regexs)
    plot_occurrences(dict, dict_2, dict_3)
    path_fastp = argv[2]
    diff = subprocessing(path_fastp, fastq)
    if not diff:
        print("error accrued")
    else:
        print(diff)
if __name__ == '__main__':
    if len(argv) == 3:
        main(argv)