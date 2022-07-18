import numpy as np
from NCellShell import *

S_i = []
SD_i = []
D_i = []
N_l = []
for i in range(9):
    S_i.append([0])
    SD_i.append([0])
    D_i.append([0])
    N_l.append(500)
target_alts = [500]
alt_edges = np.linspace(300, 600, num=10)
print(alt_edges)
lam = [1000]

my_cell = NCell(S_i, SD_i, D_i, N_l, target_alts, alt_edges, lam)
my_cell.save("./", "test_save")
cell_cpy = NCell.load("./test_save/")
print("Success")