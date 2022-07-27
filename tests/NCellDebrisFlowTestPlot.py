import sys
sys.path.append('./../')

from NCellShell import NCell
 
atmosphere = NCell.load("./data/NCellDebrisFlowTest/")
t = atmosphere.get_t()
T = t[-1]
N = atmosphere.get_N()

import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()
ax1.set_xlabel('time (yr)')
ax1.set_ylabel('log(number)')
ax1.set_yscale('log')
for i in range(0, len(N)):
    ax1.plot(t, N[i], label='SD'+str(i))
ax1.set_ylim(1, 1e5)
ax1.set_xlim(0,T)
ax1.legend()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
