from utils import GreedyMotifSearch, RandomizedMotifSearch, GibbsSampler, Consensus, Score


if __name__ == '__main__':
    DosR = []
    with open('data//DosR.txt') as file:
        DosR = file.read().splitlines()
    print(DosR)
    motifs = GreedyMotifSearch(DosR, 15, WithPseudocounts=False)
    consensus = Consensus(motifs)
    print(consensus)
    print(Score(motifs), consensus)
    motifs = RandomizedMotifSearch(DosR, 15)
    consensus = Consensus(motifs)
    print(consensus)
    print(Score(motifs), consensus)
    motifs = GibbsSampler(DosR, 15, 10, 100)
    consensus = Consensus(motifs)
    print(motifs)
    print(consensus)
    print(Score(motifs), consensus)
