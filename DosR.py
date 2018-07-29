from utils import GreedyMotifSearch, Consensus


if __name__ == '__main__':
    DosR = []
    with open('data//DosR.txt') as file:
        DosR = file.read().splitlines()
    print(DosR)
    motifs = GreedyMotifSearch(DosR, 15)
    print(motifs)
    print(Consensus(motifs))
