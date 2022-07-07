
from NCellShell import *

S_i = [[0]]
SD_i = [[0]]
D_i = [[0]]
N_l = [500]
target_alts = [500]
alt_edges = [475, 525]
lam = [1000]

my_cell = NCell(S_i, SD_i, D_i, N_l, target_alts, alt_edges, lam)
my_cell.save("./", "test_save")
cell_cpy = NCell.load("./test_save/")
print("Success")