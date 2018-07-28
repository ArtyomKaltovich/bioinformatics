from utils import MinimumSkew

if __name__ == '__main__':
    with open('data//e_coli.txt') as file:
        e_coli = file.read()
    print(MinimumSkew(e_coli))