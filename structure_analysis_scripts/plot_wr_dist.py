import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('ticks')

lks = range(0,4)

for lk in lks:
    wr = np.loadtxt('data/LK_%s/wr' %lk)
    sns.distplot(wr)

plt.show()
