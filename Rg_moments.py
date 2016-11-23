#Calculate the  principle moments of the radius of gyration of a polymer represented by 
#a vector of positions
from numpy import *
def Rg_moments(r):
    N=len(r[:,0])
    #calculate the center of mass
    Rcomx=sum(r[:,0])/N
    Rcomy=sum(r[:,1])/N
    Rcomz=sum(r[:,2])/N
    Rgxx=sum((r[:,0]-Rcomx)*(r[:,0]-Rcomx))
    Rgyy=sum((r[:,1]-Rcomy)*(r[:,1]-Rcomy))
    Rgzz=sum((r[:,2]-Rcomz)*(r[:,2]-Rcomz))
    Rgxy=sum((r[:,0]-Rcomx)*(r[:,1]-Rcomy))
    Rgxz=sum((r[:,0]-Rcomx)*(r[:,2]-Rcomz))
    Rgyz=sum((r[:,1]-Rcomy)*(r[:,2]-Rcomz))
    Rgyx=Rgxy
    Rgzy=Rgyz
    Rgzx=Rgxz
    Rg=array([[Rgxx,Rgxy,Rgxz],[Rgyx,Rgyy,Rgyz],[Rgzx,Rgzy,Rgzz]])/N
    Rg_moments=linalg.eigvals(Rg)
    Rg_moments=sort(Rg_moments)
    return Rg_moments

    
