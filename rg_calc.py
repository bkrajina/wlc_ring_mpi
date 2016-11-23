from numpy import *
from Rg_moments_stat import *

folder='data'
nmax=1000

moments=Rg_moments_stat(folder,nmax)[0,:]

rg=[sum(moments)]
rgrms=[sum(moments)**0.5]
savetxt('data/RGYSQ_AVG',rg)
savetxt('data/RGRMS_AVG',rgrms)
