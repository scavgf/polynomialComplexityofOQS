include("./parameter.jl")
using LinearAlgebra
using FileIO, JLD2
using .Parameters:Omega,alphac,omegac, sigmax,sigmay,sigmaz,sigma0,Dt,Nt,lambda,dim,K,epsilon

##initialization
rho0=[1/2,0,0,1/2]
Sxn=zeros(Float64,4,4)
Sxn[3,4]=-2
Sxn[4,3]=2
#free evolution 
##system's Hamiltonian Hs=Omega*sigmax
Hsn=Omega.*Sxn

file1=open("data/Szt-D-"*string(dim)*"-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-lambda="*string(lambda)*"-HO.dat","w")
file2=open("data/BDt-D-"*string(dim)*"-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-lambda="*string(lambda)*"-HO.dat","w")
for ix=10:10:Nt
        namet="tmp/GT-"*string(ix)*"-dim="*string(dim)*"-alpha="*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*".jld"
        if isfile(namet)
            data=load(namet)
            Gt=data["Gt"]
            dimc=data["dimc"]
            t=ix*Dt
            rhot=exp(t*Hsn)*Gt*rho0
		#################################
		    Sz=[0,0,0,1]'*rhot
		    S0=[1,0,0,0]'*rhot
            write(file1, join([t, Sz], " "), "\n")
            write(file2, join([t, dimc], " "), "\n")
		    println("ix=", ix, " dimc=", dimc, " Sz=", Sz, " error=",1-2*S0)
        else 
            break
        end
end
