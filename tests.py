import unittest
from time import process_time

from utils import *


class utils_test(unittest.TestCase):
    def test_PatternMatching(self):
        c = PatternMatching("ATAT","GATATATGCATATACTT")
        correct = [1, 3, 9]
        self.assertEqual(correct, c)
        c = PatternMatching("CTTGATCAT","CTTGATCATCTTGATCATCTTGATCAT")
        correct = [0, 9, 18]
        self.assertListEqual(correct, c)


    def test_FrequencyMap(self):
        c = FrequencyMap("CGATATATCCATAG", 3)
        correct = {'CGA': 1, 'GAT': 1, 'ATA': 3, 'TAT': 2, 'ATC': 1, 'TCC': 1, 'CCA': 1, 'CAT': 1, 'TAG': 1}
        self.assertDictEqual(c, correct)

    def test_FrequentWords(self):
        c = FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
        correct = ["CATG", "GCAT"]
        self.assertDictEqual(Counter(c), Counter(correct))  # compare lists order free O(n)

    def test_Reverse(self):
        c = Reverse("AAAACCCGGT")
        correct = "TGGCCCAAAA"
        self.assertEqual(correct, c)

    def test_Complement(self):
        c = Complement("AAAACCCGGT")
        correct = "TTTTGGGCCA"
        self.assertEqual(correct, c)

    def test_ReverseComplement(self):
        c = ReverseComplement("AAAACCCGGT")
        correct = "ACCGGGTTTT"
        self.assertEqual(correct, c)

    def test_fast_and_slow_SymbolArray(self):
        c = SymbolArray("AAAAGGGG", "A")
        correct = {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
        self.assertDictEqual(correct, c)
        c2 = FasterSymbolArray("AAAAGGGG", "A")
        self.assertDictEqual(c, c2)
        c = FasterSymbolArray("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "CC")
        correct = {0: 7, 1: 7, 2: 7, 3: 7, 4: 7, 5: 7, 6: 7, 7: 7, 8: 7, 9: 7, 10: 7, 11: 7, 12: 7, 13: 7, 14: 7, 15: 7, 16: 7, 17: 7, 18: 7, 19: 7, 20: 7, 21: 7, 22: 7, 23: 7, 24: 7, 25: 7, 26: 7, 27: 7, 28: 7, 29: 7, 30: 7, 31: 7, 32: 7, 33: 7, 34: 7, 35: 7, 36: 7, 37: 7, 38: 7, 39: 7, 40: 7, 41: 7, 42: 7, 43: 7, 44: 7, 45: 7, 46: 7, 47: 7, 48: 7, 49: 7, 50: 7, 51: 7, 52: 7, 53: 7, 54: 7, 55: 7, 56: 7, 57: 7, 58: 7, 59: 7, 60: 7, 61: 7, 62: 7, 63: 7, 64: 7, 65: 7, 66: 7, 67: 7, 68: 7, 69: 7, 70: 7, 71: 7, 72: 7, 73: 7, 74: 7, 75: 7, 76: 7, 77: 7, 78: 7, 79: 7, 80: 7, 81: 7, 82: 7, 83: 7, 84: 7, 85: 7, 86: 7, 87: 7, 88: 7, 89: 7, 90: 7, 91: 7, 92: 7, 93: 7, 94: 7, 95: 7, 96: 7, 97: 7, 98: 7, 99: 7, 100: 7, 101: 7, 102: 7, 103: 7, 104: 7, 105: 7, 106: 7, 107: 7, 108: 7, 109: 7, 110: 7, 111: 7, 112: 7, 113: 7, 114: 7, 115: 7, 116: 7, 117: 7, 118: 7, 119: 7, 120: 7, 121: 7, 122: 7, 123: 7, 124: 7, 125: 7, 126: 7, 127: 7, 128: 7, 129: 7, 130: 7, 131: 7, 132: 7, 133: 7, 134: 7}
        self.assertDictEqual(correct, c)

    def test_SkewArray(self):
        c = SkewArray("CATGGGCATCGGCCATACGCC")
        correct = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
        self.assertListEqual(correct, c)

    def test_SkewArray_raise_error_on_invalid_param(self):
        with self.assertRaises(ValueError):
            SkewArray("123")

    def test_MinimumSkew(self):
        correct = [11, 24]
        c = MinimumSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
        self.assertListEqual(correct, c)

    def test_HammingDistance(self):
        correct = 3
        c = HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC")
        self.assertEqual(correct, c)
        c = HammingDistance("123", "123")
        self.assertEqual(0, c)
        c = HammingDistance("qwe", "123")
        self.assertEqual(3, c)

    def test_HammingDistance_on_empty_strings(self):
        c = HammingDistance("", "")
        self.assertEqual(0, c)

    def test_HammingDistance_on_different_length_string(self):
        c = HammingDistance("123", "2dddddd")
        self.assertEqual(7, c)
        c = HammingDistance("123456", "qw3")
        self.assertEqual(5, c)

    def test_ApproximatePatternMatching(self):
        correct = [6, 7, 26, 27]
        c = ApproximatePatternMatching("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT",
                                       "ATTCTGGA", 3)
        self.assertEqual(correct, c)

    def test_ApproximatePatternCount(self):
        c = ApproximatePatternCount("GAGG", "TTTAGAGCCTTCAGAGG", 2)
        self.assertEqual(4, c)

    def test_CountMotifs(self):
        correct = {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}
        c = CountMotifs(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"])
        self.assertDictEqual(correct, c)

    def test_ProfileMotifs(self):
        correct = {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4], 'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0], 'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2], 'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}
        c = ProfileMotifs(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"])
        self.assertDictEqual(correct, c)

    def test_Consensus(self):
        correct = "CACCTA"
        c = Consensus(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"])
        self.assertEqual(correct, c)

    def test_Score(self):
        correct = 14
        c = Score(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"])
        self.assertEqual(correct, c)

    def test_Score_precomputed_consensus(self):
        correct = 14
        c = Score(["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], "CACCTA")
        self.assertEqual(correct, c)

    def test_Pr(self):
        correct = 0.0008398080000000002
        profile = { "A": [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
                    "C": [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
                    "G": [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
                    "T": [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}
        c = Pr("ACGGGGATTACC", profile)
        self.assertEqual(correct, c)

    def test_ProfileMostProbableKmer(self):
        correct = "CCGAG"
        profile = { "A": [0.2, 0.2, 0.3, 0.2, 0.3],
              "C": [0.4, 0.3, 0.1, 0.5, 0.1],
              "G": [0.3, 0.3, 0.5, 0.2, 0.4],
              "T": [0.1, 0.2, 0.1, 0.1, 0.2]
            }
        c = ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, profile)
        self.assertEqual(correct, c)
        correct = "CAGCG"
        profile = { "A": [0.2, 0.2, 0.3, 0.2, 0.3],
                    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
                    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
                    "T": [0.1, 0.2, 0.1, 0.1, 0.2]}
        c = ProfileMostProbableKmer("TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA", 5, profile)
        self.assertEqual(correct, c)

    def test_ProfileMostProbableKmer_ties(self):
        correct = "AAC"
        profile = { "A": [0.1, 0.1, 0.1],
              "C": [0.0, 0.0, 0.0],
              "G": [0.0, 0.0, 0.0],
              "T": [0.0, 0.0, 0.0]
            }
        c = ProfileMostProbableKmer("AACCGGTT", 3, profile)
        self.assertEqual(correct, c)

    def test_GreedyMotifSearch(self):
        correct = ["CAG", "CAG", "CAA", "CAA", "CAA"]
        dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
        c = GreedyMotifSearch(dna, 3)
        self.assertListEqual(correct, c)

    def test_Entropy(self):
        correct = 0.971
        profile = {'A': [0.0], 'C': [0.6], 'G': [0.0], 'T': [0.4]}
        c = Entropy(profile)
        self.assertAlmostEqual(correct, c, delta=0.01)

        correct = 9.916290005356972
        profile = {
            'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
            'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
            'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
            'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
        }
        c = Entropy(profile)
        self.assertAlmostEqual(correct, c, delta=0.01)
