from collections import Counter
import math


ComplementMap = {"A": "T", "C": "G", "G": "C", "T": "A"}


def PatternCount(Pattern, Text: str):
    return Text.count(Pattern)


def PatternMatching(Pattern, Genome):
    """
    Pattern Matching Problem:  Find indexes of all occurrences of a pattern in a Genome.
    :return: All starting positions in Genome where Pattern appears as a substring.
    :examples: PatternMatching("ATAT","GATATATGCATATACTT") -> 1 3 9
    """
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions


def frequency_map(Text, k):
    """
    generates a frequency map from a given string Text and integer k.

    :examples:
        For Text = "CGATATATCCATAG" result will be:
        ATA --> 3
        ATC --> 1
        CAT --> 1
        CCA --> 1
        CGA --> 1
        GAT --> 1
        TAT --> 2
        TCC --> 1
        TAG --> 1
    """
    c = Counter()
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        c[Pattern] += 1
    return c


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq


def FrequentWords(Text, k):
    freq = frequency_map(Text, k)
    m = max(freq.values())
    words = [key for key in freq if freq[key] == m]
    return words


def Reverse(Pattern):
    result = ''
    for i in range(len(Pattern) - 1, -1, -1):
        result += Pattern[i]
    return result


def Complement(Pattern):
    result = map(lambda x: ComplementMap[x], Pattern)
    return "".join(result)


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern


def SymbolArray(Genome, symbol):
    """
    We will keep track of the total number of occurrences of C that we encounter in each window of ExtendedGenome
        by using a symbol array. The i-th element of the symbol array is equal to the number of occurrences
        of the symbol in the window of length len(Genome)//2 starting at position i of ExtendedGenome.
    :param Genome:
    :param symbol:
    :return:
    :example:
        Sample Input:
            AAAAGGGG
            A
        Sample Output:
            {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def SkewArray(Genome: str):
    """
    We will keep track of the difference between the total number of occurrences of G and the total number of
        occurrences of C that we have encountered so far in Genome by using a skew array. This array, denoted Skew,
        is defined by setting Skew[i] equal to the number of occurrences of G minus the number of occurrences of C
        in the first i nucleotides of Genome. We also set Skew[0] equal to zero.

        Given a string Genome, we can form its skew array by setting Skew[0] equal to 0,
        and then ranging﻿ through the genome.  At position i of Genome, if we encounter an A or a T,
        we set Skew[i+1] equal to Skew[i]; if we encounter a G, we set Skew[i+1] equal to Skew[i]+1;
        if we encounter a C, we set Skew[i+1] equal to Skew[i]-1.
    :return: SkewArray (see function description)
    :example:
        Sample Input:
            CATGGGCATCGGCCATACGCC
        Sample Output:
            0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    """

    genome_length = len(Genome)
    skew = [0,] * (genome_length + 1)
    for i in range(genome_length):
        if Genome[i] in ['g', 'G']:
            skew[i + 1] = skew[i] + 1
        elif Genome[i] in ['c', 'C']:
            skew[i + 1] = skew[i] - 1
        elif Genome[i] in ['t', 'T', 'a', 'A']:
            skew[i + 1] = skew[i]
        else:
            raise ValueError("you should pass correct genome here (only letters a,c,g,t in any case are allowed)")
    return skew


def MinimumSkew(Genome):
    """
    generate an empty list positions
        set a variable equal to SkewArray(Genome)
        find the minimum value of all values in the skew array
        range over the length of the skew array and add all positions achieving the min to positions
    :param Genome: A DNA string. 
    :return: All integer(s) i minimizing Skew[i] among all values of i (from 0 to len(Genome)).
    :example:
        Sample Input:
            TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT
        Sample Output:
            11 24 (indexes of all minimal items)
    """
    result = SkewArray(Genome)
    min_item = result[0]
    min_index = [0]
    for index, x in enumerate(result[1:]):
        if x < min_item:
            min_index = [index + 1]
            min_item = x
        elif x == min_item:
            min_index.append(index + 1)
    return min_index


def HammingDistance(p, q):
    """
    Compute the Hamming distance between two strings.
    :return: The Hamming distance between these strings.
    :example:
        Sample Input:
            GGGCCGTTGGT
            GGACCGTTGAC
        Sample Output:
            3
    """
    result = 0
    for x, y in zip(p, q):
        if x != y:
            result += 1
    return result + abs(len(p) - len(q))


def ApproximatePatternMatching(Genome, Pattern, d=3):
    """
    Find all approximate occurrences of a pattern in a string.
        We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some
        k-mer substring Pattern' of Text having d or fewer mismatches with Pattern; that is,
        HammingDistance(Pattern, Pattern') ≤ d.
        Our observation that a DnaA box may appear with slight variations leads to the following generalization
        of the Pattern Matching Problem.
    :param Genome: Genome in which search is occur.
    :param Pattern: Pattern which is been looking for.
    :param d: Maximal Hamming distance between string to call them approximately equal.
    :return:  All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    :example:
        Sample Input:
            CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
            ATTCTGGA
            3
        Sample Output:
            6 7 26 27
    """
    positions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if HammingDistance(Genome[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Pattern, Genome, d=3):
    """
    Given input strings Text and Pattern as well as an integer d, we extend the definition of PatternCount
        to the function ApproximatePatternCount(Pattern, Text, d).
        This function computes the number of occurrences of Pattern in Text with at most d mismatches.
    :example:
        Sample Input:
            GAGG
            TTTAGAGCCTTCAGAGG
            2
            Note! Param order is different from ApproximatePatternMatching function due to course source checking system
        Sample Output:
            4
    """
    count = 0
    for i in range(len(Genome) - len(Pattern) + 1):
        if HammingDistance(Genome[i:i+len(Pattern)], Pattern) <= d:
            count = count + 1
    return count


def CountMotifs(motifs: list) -> dict:
    """
    Function takes a list of strings Motifs as input and returns the count matrix of Motifs (as a dictionary of lists).
    :param motif: a list of strings which represent a genome part
    :return: dictionary where keys are A, C, G, T and values are list with their occurrences in patterns
        on that index.
    :example:
    Sample Input:
        AACGTA
        CCCGTT
        CACCTT
        GGATTA
        TTCCGG
    Sample Output:
        {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}
    """
    result = {letter: [0] * len(motifs[0]) for letter in "ACGT"}
    for i in range(len(motifs[0])):
        for motif in motifs:
            result[motif[i]][i] += 1
    return result


def ProfileMotifs(Motifs: list) -> dict:
    """
    Function that takes Motifs as input and returns their profile matrix as a dictionary of lists.
    :param Motifs: A list of kmers Motifs
    :return: the profile matrix of Motifs, as a dictionary of lists.
    :example:
    Sample Input:
        AACGTA
        CCCGTT
        CACCTT
        GGATTA
        TTCCGG
    Sample Output:
        {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4], 'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0], 'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2], 'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}
    """
    result = {key: [i / len(Motifs) for i in value] for key, value in CountMotifs(Motifs).items()}
    return result


def Consensus(Motifs: list) -> str:
    """
    Form a consensus string, from the most popular nucleotides in each column of the motif matrix
        (ties are broken arbitrarily). If we select Motifs correctly from the collection of upstream regions,
        then Consensus(Motifs) provides a candidate regulatory motif for these regions.
    :param Motifs: A set of kmers Motifs
    :return: A consensus string of Motifs.
    :example:
    Sample Input:
        AACGTA
        CCCGTT
        CACCTT
        GGATTA
        TTCCGG
    Sample Output:
        CACCTA
    """
    k = len(Motifs[0])
    count = CountMotifs(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs: list, consensus: str=None) -> list:
    """
    Constructing Consensus(Motifs) and then summing the number of symbols in the j-th column of Motifs
        that do not match the symbol in position j of the consensus string.
    :param Motifs: A set of kmers Motifs
    :param Consensus: Precomputed Consensus, can be ommited
    :return: The score of these k-mers.
    :example:
    Sample Input:
        AACGTA
        CCCGTT
        CACCTT
        GGATTA
        TTCCGG
    Sample Output:
        14
    """
    consensus = Consensus(Motifs) if consensus is None else consensus
    result = 0
    for motif in Motifs:
        for i, letter in enumerate(motif):
            if letter != consensus[i]:
                result += 1
    return result


def Pr(Text: str, Profile: dict) -> float:
    """
    Calculate probability of Text by profile matrix
    :example:
    Sample Input:
        ACGGGGATTACC
        0.2 0.2 0.0 0.0 0.0 0.0 0.9 0.1 0.1 0.1 0.3 0.0
        0.1 0.6 0.0 0.0 0.0 0.0 0.0 0.4 0.1 0.2 0.4 0.6
        0.0 0.0 1.0 1.0 0.9 0.9 0.1 0.0 0.0 0.0 0.0 0.0
        0.7 0.2 0.0 0.0 0.1 0.1 0.0 0.5 0.8 0.7 0.3 0.4
    Sample Output:
        0.0008398080000000002
    """
    letter_to_row = {key: i for i, key in enumerate("ACGT")}
    result = 1
    for index, letter in enumerate(Text):
        result *= Profile[letter][index]
    return result


def ProfileMostProbableKmer(text, k, profile):
    """
    Find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile
        among all k-mers in Text. If there are multiple Profile-most probable k-mers in Text,
        then we select the first such k-mer occurring in Text.
    :example:
    Sample Input:
        ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
        5
        { A: 0.2 0.2 0.3 0.2 0.3
          C: 0.4 0.3 0.1 0.5 0.1
          G: 0.3 0.3 0.5 0.2 0.4
          T: 0.1 0.2 0.1 0.1 0.2
        }
    Sample Output:
        CCGAG
    """
    maximim = (text[0:k], Pr(text[0:k], profile))
    for i in range(1, len(text) - k + 1):
        probability = Pr(text[i:i + k], profile)
        if probability > maximim[1]:
            maximim = (text[i:i + k], probability)
    return maximim[0]


def GreedyMotifSearch(Dna, k, t=None):
    """
    :input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
    :example:
    Sample Input:
        GGCGTTCAGGCA
        AAGAATCAGTCA
        CAAGGAGTTCGC
        CACGTCAATCAC
        CAATAATATTCG
        3 5
    Sample Output:
        CAG
        CAG
        CAA
        CAA
        CAA
    """
    t = t if t is not None else len(Dna)
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileMotifs(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def Entropy(profile: dict) -> float:
    """
    Calculates the entropy of profile matrix.
        H(p1,…,pN)=−∑i=1Npi⋅log2pi
    """
    entropy_log = lambda x : x * math.log(x, 2) if x != 0 else 0  # log(0, 2) isn't defined but assumes as 0
                                                                  # in entropy calculations
    result = [[entropy_log(x) for row in list(profile.values()) for x in row ]]
    result = sum([sum(row) for row in result])
    return -result
