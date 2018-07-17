from collections import Counter


ComplementMap = {"A": "T", "C": "G" , "G": "C", "T": "A"}


def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


def PatternMatching(Pattern, Genome):
    """
    Pattern Matching Problem:  Find all occurrences of a pattern in a Genome.
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
