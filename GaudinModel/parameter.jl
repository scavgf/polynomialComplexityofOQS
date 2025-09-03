module Parameters
using LinearAlgebra
#parameters for zero temperature case
#Pauli matrix 
###########################################################
sigmax=[0 1; 1 0]
sigmay=[0 -im; im 0]
sigmaz=[1 0; 0 -1]
sigma0=[1 0; 0 1]
#Generation and aninihilation opeator 
###########################################################
#super-operators
I0=[1.0 0 0 0; 0 1.0 0 0; 0 0 1.0 0; 0 0 0 1.0]
###########################################################
Sxp=[0 1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0]
Syp=[0 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0]
Szp=[0 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0]
Sxn=[0 0 0 0; 0 0 0 0; 0 0 0 -2; 0 0 2 0]
Syn=[0 0 0 0; 0 0 0 2; 0 0 0 0; 0 -2 0 0]
Szn=[0 0 0 0; 0 0 -2 0; 0 2 0 0; 0 0 0 0]
###########################################################
Sxxpp=[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]
Syypp=[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0]
Szzpp=[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]
Sxypp=[0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0]
Syzpp=[0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0]
Szxpp=[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 1 0 0]
###########################################################
Sxypn=[0 0 0 2; 0 0 0 0; 0 0 0 0; 0 0 0 0]
Syzpn=[0 2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
Szxpn=[0 0 2 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
###########################################################
const Omega=0.0/2 #energy splitting 
const A=1.0/4.0 #the interaction strength
const N=49 #number of nucler spin
const Od=2 #order of method
const Dt=0.03
const Nt=500 #evoution steps
const dimG=7 #dimension of generator
const dims=4 #dimension of spin
const dim=100 #bond-dimension
const lambda=1.0 #bond-dimension
const epsilon=1.0e-6
const rhoB0=[1/2,0,0,0]
const Tr=[2 0 0 0]
###########################################################
end
