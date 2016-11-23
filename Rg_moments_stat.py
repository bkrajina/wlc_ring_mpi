#Calculate the average, standard deviation, and standard error of the principle moment of 
#a collection of saved polymer configurations contained within a root data folder
#Configuration files are assumed to have a labelling scheme like 'rn' where n is an integer
#Save the statistics in the root folder
#The order determines whether to calculate averages and standard errors  based on Rg^2 or sqrt(Rg^2)
#Note that if order=1 is selected, only the standard error of the mean is meaningful, since we are computing <Rg^2>^0.5 not <Rg>

from numpy import *
from Rg_moments import *


def Rg_moments_stat(folder,nmax,n0=1,nskip=1,order=2):
    index=n0
    i=0
    #Number of total saved configurations
    nsave=1+(nmax-n0)/nskip
    Rg_moments_all=zeros([nsave,3])
    while (index<=nmax):
        r=loadtxt(folder+'/r'+str(index))
        Rg_moments_all[i,:]=Rg_moments(r)
        i=i+1
        index=index+nskip
    Rg1avg=average(Rg_moments_all[:,0])
    Rg2avg=average(Rg_moments_all[:,1])
    Rg3avg=average(Rg_moments_all[:,2])
    Rg1stdv=std(Rg_moments_all[:,0])
    Rg2stdv=std(Rg_moments_all[:,1])
    Rg3stdv=std(Rg_moments_all[:,2])
    Rg1stdv=Rg1stdv/(2*(Rg1avg**0.5))
    Rg2stdv=Rg2stdv/(2*(Rg2avg**0.5))
    Rg3stdv=Rg3stdv/(2*(Rg3avg**0.5))
    Rg1stderr=Rg1stdv/sqrt(nsave)
    Rg2stderr=Rg2stdv/sqrt(nsave)
    Rg3stderr=Rg3stdv/sqrt(nsave)
    if (order==1):
        Rg1avg=Rg1avg**0.5
        Rg2avg=Rg2avg**0.5
        Rg3avg=Rg3avg**0.5

        Rg1stderr=Rg1stderr/(2*(Rg1stderr**0.5))
        Rg2stderr=Rg2stderr/(2*(Rg2stderr**0.5))
        Rg3stderr=Rg3stderr/(2*(Rg3stderr**0.5))
 
    data=zeros([3,3])
    data[0,:]=array([Rg1avg,Rg2avg,Rg3avg])
    data[1,:]=array([Rg1stdv,Rg2stdv,Rg3stdv])
    data[2,:]=array([Rg1stderr,Rg2stderr,Rg3stderr])
    return data
