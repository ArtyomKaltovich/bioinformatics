import matplotlib.pyplot as plt
from utils import FasterSymbolArray

with open('data//e_coli.txt') as file:
    e_coli = file.read()

array = FasterSymbolArray(e_coli, "C")


plt.plot(*zip(*sorted(array.items())))
plt.title("Amount of Cytosine in the Following Nucleotides")
plt.xlabel("Letter Number")
plt.ylabel("#Cytosine")
plt.savefig('plots//e_coli_c.png')