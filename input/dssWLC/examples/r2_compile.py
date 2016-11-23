from numpy import *

L=[0.125,0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0]
R2=[]
for l in L:
    R2.append(loadtxt('R2_avg'+'_L_' + str(l)))
data=zeros([len(L),2])
data[:,0]=L
data[:,1]=R2
savetxt('r2_vs_L_lena',data)
