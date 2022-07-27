import sys
sys.path.append('./../')

from NCellShell import NCell
 
atmosphere = NCell.load("./data/NCellDerelictFlowTest/")
t = atmosphere.get_t()
T = t[-1]
D = atmosphere.get_D()

import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()
ax1.set_xlabel('time (yr)')
ax1.set_ylabel('log(number)')
ax1.set_yscale('log')
for i in range(0, len(D)):
    ax1.plot(t, D[i][0], label='D'+str(i))
ax1.set_ylim(1, 1e2)
ax1.set_xlim(0,T)
ax1.legend()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
