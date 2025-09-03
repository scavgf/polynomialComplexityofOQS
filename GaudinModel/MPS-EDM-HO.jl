include("./parameter.jl")
include("./EV-MPO-HO-canonical.jl")
using LinearAlgebra
using TensorOperations 
using DelimitedFiles
using ZChop
using Dates
using FileIO, JLD2
using Polynomials
using .Evolution:evolution, propagator
using .Parameters:Sxp,Syp,Szp,Sxn,Syn,Szn,I0,dimG
using .Parameters:Omega, Dt,Nt,dim,N,A,epsilon,Od,lambda
#super operator in terms of [I, sx,sy,sz]
########################
#free evolution 
##system's Hamiltonian Hs=Omega*sigmax
Hsn=Omega.*Szn
#################configruation#################
#  B-B---B-B-----B-B----B-B (correlation tensor)
#  | |   | |     | |    | |
#  G-G----G-G----G-G----G-G  (progator tenor)
# generating MPS
###########################
function MPSL(ix::Int64, coe::ComplexF64)
	G=zeros(ComplexF64, 4, dimG, 4)
	#interaction picture
	t=(ix+0.5)*Dt
	Sxnt=exp(-t*Hsn)*Sxn*exp(t*Hsn)
	Sxpt=exp(-t*Hsn)*Sxp*exp(t*Hsn)
	Synt=exp(-t*Hsn)*Syn*exp(t*Hsn)
	Sypt=exp(-t*Hsn)*Syp*exp(t*Hsn)
	#############################
	G[:,1,:]=I0
	G[:,2,:]=Dt*lambda*coe.*Sxpt
	G[:,3,:]=Dt*lambda*coe.*Sypt
	G[:,4,:]=Dt*lambda*coe.*Szp
	G[:,5,:]=Dt*lambda*coe.*Sxnt
	G[:,6,:]=Dt*lambda*coe.*Synt
	G[:,7,:]=Dt*lambda*coe.*Szn
	return G
end
#########################################
println(" Dt=", Dt," Nt=", Nt, " dim=",dim, " epsilon=", epsilon)
##############generating MPS########################
ply=Polynomial([1.0/factorial(i) for i=0:Od])
println("ply=",ply)
Roots=-1.0./roots(ply)
println("roots=",Roots)
##############################################
MPS=Vector{Array{ComplexF64,3}}()
for ix=Nt:-1:1
	for coe in Roots
		L=MPSL(ix,coe)
		append!(MPS,[L])
	end
end
##############################################
##picking tensor
P=zeros(dimG,dimG,dimG)
P[1,:,:]=Matrix{Float64}(I,dimG,dimG)
for ix=2:dimG
	P[ix,1,ix]=1.0
end
#############B tensor###########################
function MPOT(Jk::Float64)
	B=zeros(4,dimG,4)
	B[:,1,:]=I0
	B[:,2,:]=Jk/lambda.*Sxn
	B[:,3,:]=Jk/lambda.*Syn
	B[:,4,:]=Jk/lambda.*Szn
	B[:,5,:]=Jk/lambda.*Sxp
	B[:,6,:]=Jk/lambda.*Syp
	B[:,7,:]=Jk/lambda.*Szp
	@tensor T[l,u,d,r]:=B[l,x,r]*P[x,u,d]
	return T
end
####BTensor###################################
println("time now=", now())
for k=1:N #N bath spin
	global MPS
	Nx=(N+1-k)
	Jk=A*sqrt(6*N/(2*N^2+3*N+1))*Nx/N
	println("Jk=",Jk)
	Tk=MPOT(Jk)
	MPS,BD=evolution(MPS, Tk)
	println("k=",k," maximal dimension=", last(BD,30))
	##############################################
	Lt=propagator(MPS, Roots[Od])
	file=open("data/propagator-Omega="*string(2*Omega)*"-N="*string(k)*"-D-"*string(dim)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-lambda="*string(lambda)*"-ep="*string(epsilon)*"-Od="*string(Od)*".dat","w")
	for iy=1:Nt-1
		write(file, join(real(Lt[iy]), " "), "\n")
	end
	close(file)
	###############################################
	fileBD=open("data/BD-Omega="*string(2*Omega)*"-N="*string(k)*"-D-"*string(dim)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-lambda="*string(lambda)*"-ep="*string(epsilon)*"-Od="*string(Od)*".dat","w")
	write(fileBD, join(BD, " "), "\n")
	close(fileBD)
	##############################################
	println("time now=", now())
	println("****************************")
end
##########save data###########################
