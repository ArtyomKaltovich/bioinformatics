import unittest
from collections import Counter

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