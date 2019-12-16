from numpy import *

n0=1
nmax=1000
dr=zeros([nmax])
L=0.5
for i in range(1,nmax+1):
    r=loadtxt('data/r' + str(i))
    dr[i-1]=sum((r[-1,:]-r[0,:])**2.)**0.5
dr_avg=[average(dr)]
dr_norm=[average(dr)/L]
savetxt('data/dr_avg',dr_avg)
savetxt('data/dr_norm',dr_norm)
