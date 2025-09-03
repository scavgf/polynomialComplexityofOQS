module Parameters
using LinearAlgebra
#parameters for zero temperature case
Omega=0.5 #energy splitting 
omegac=5.0
#Pauli matrix 
sigmax=[0 1; 1 0]
sigmay=[0 -im; im 0]
sigmaz=[1 0; 0 -1]
sigma0=[1 0; 0 1]
#Generation and aninihilation opeator 
###########################################################
#const Dt=0.01#evolution step
#const Nt=400#evolutin time Nt*Dt
#const K=500
#const dim=80 #bond-dimension truncation
#const lambda=1.0 #measurement strength
#const epsilon=1.0e-4#measurement strength
alphac=parse(Float64, ARGS[1])
Dt=parse(Float64, ARGS[2])
Nt=parse(Int, ARGS[3])
K=parse(Int, ARGS[4])
dim=parse(Int, ARGS[5])
epsilon=parse(Float64, ARGS[6])
lambda=parse(Float64, ARGS[7])
###########################################################
end
