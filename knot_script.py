from knot_check import *
from numpy import *
from scipy import *
fileroot='data/r'
savefolder='data/'
Nmin=1;
Nmax=10;
delta=zeros([Nmax-Nmin+1,1])

for fileindex in range(Nmin,Nmax+1):
    delta[fileindex-Nmin,0]=alexander(fileroot+str(fileindex),savefolder[0:-1])

delta=rint(delta)
delta=delta.astype(int)
savetxt(savefolder + 'delta',delta.astype(str),fmt="%s")

