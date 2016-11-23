import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc

rc('axes',labelsize=22.)
rc('xtick',labelsize=18.)
rc('ytick',labelsize=18.)

sns.set_style('ticks')
Lks = range(0,6)
lk0 = (1.e4)/10.5
Rgs = []


for lk in Lks:
    Rgs.append(np.loadtxt('data/LK_%s/RGYRSQ_Avg' %(lk))**0.5)

sig = np.array(Lks)/lk0

plt.plot(sig,Rgs,'-ro')

plt.show()
