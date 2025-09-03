include("./parameter.jl")
include("./EV-MPO-GT-HO.jl")
using LinearAlgebra
using TensorOperations 
using Polynomials
using DelimitedFiles
using ZChop
using FileIO, JLD2
using .Evolution:evolution, propagator, bondDimension
using .Parameters:Omega,alphac,omegac, sigmax,sigmay,sigmaz,sigma0,Dt,Nt,lambda,dim,K,epsilon
#super operator in terms of [I, sx,sy,sz]
Sxn=zeros(Float64,4,4)
Sxn[3,4]=-2
Sxn[4,3]=2
Szn=zeros(Float64,4,4)
Szn[2,3]=-2
Szn[3,2]=2
Szp=zeros(Float64,4,4)
Szp[1,4]=1
Szp[4,1]=1
#free evolution 
##system's Hamiltonian Hs=Omega*sigmax
Hsn=Omega.*Sxn
##decoherence induced  by eqal-time correlation
Dtx=omegac*Dt
#generator 
dimG=3 #dimension of generator
dimB=2  #bond dimension of correlation tensor
#############MPS tensor #################
println("alpha=",alphac," Dt=", Dt," Nt=", Nt, " K=",K," dim=",dim, " epsilon=", epsilon)
####picking tensor#################################
IG=Matrix{Float64}(I,(dimG,dimG))
P=zeros(Float64,(dimG,dimG,dimG))
P[1,:,:]=IG
P[2,1,2]=1.0
P[3,1,3]=1.0
##############################################a
file=open("data/propagator-D-"*string(dim)*"-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-lambda="*string(lambda)*"-HO.dat","w")
fileE=open("data/entanglement-D-"*string(dim)*"-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-lambda="*string(lambda)*"-HO.dat","w")
fileD=open("data/BD-D-"*string(dim)*"-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-lambda="*string(lambda)*"-HO.dat","w")
rho0=[1/2,0,0,1/2]
############################
# generating MPS
###########################
function MPSM(ix::Int64,coe::ComplexF64)
	G=zeros(ComplexF64,(4,dimG,4))
	G[:,1,:]=Matrix{ComplexF64}(I,4,4)
	#interaction picture
	t=(ix+0.5)*Dt
	Sznt=exp(-t*Hsn)*Szn*exp(t*Hsn)
	Szpt=exp(-t*Hsn)*Szp*exp(t*Hsn)
	G[:,2,:]=Dtx*lambda*coe.*Sznt
	G[:,3,:]=Dtx*lambda*coe.*Szpt
	return G
end
############################
# generating MPO
###########################
function MPOM(ix::Int64)
	tx=ix*Dtx
	ucdt=alphac*(log(1+2*(tx^2+Dtx^2)+(tx^2-Dtx^2)^2)-2*log(1+tx^2))/(2*Dtx)^2
	usdt=alphac*(2*atan(tx)-atan(tx-Dtx)-atan(tx+Dtx))/Dtx^2
	B=zeros(Float64,(dimB,dimG,dimB))
	B[:,1,:]=Matrix{Float64}(I,(dimB,dimB))
	B[2,2,1]=ucdt
	B[2,3,1]=usdt
	@tensor T[l,u,d,r]:=P[x,u,d]*B[l,x,r]
	return T
end
####################################################
 #the expansion coefficient of high order expansion
Ol=parse(Int64, ARGS[8]) #the expansion order
ply=Polynomial([1.0/factorial(i) for i=0:Ol])
println("ply=",ply)
Roots=-1.0./roots(ply)
println("roots=",Roots)
####################################################
#generate MPS
MPS=[]
#############generate MPO T=a (b, c) d
##Tr=Picking tensor * B tensor
MPO=[]
Bl=Matrix{Float64}(I,(dimG,dimB))
Bl[1,:]=[1,0]
Bl[2,:]=[0,1]
Bl[3,:]=[0,0]
@tensor Tl[u,d,r]:=P[x,u,d]*Bl[x,r]
#boundary tensor #########
T0=MPOM(0)
MPO0=[Tl,T0,T0,T0,T0,T0]
####################################################
sl=0
for ix=0:Nt-1
	global MPS
	global MPO
	global sl
	global dimc
	#################################
	for iy=1:Ol #the expansion order to epsilon=omegac*Dt
		#println("the direction of MPS-MPS:", sl, "ix=",ix)
		###generate new MPS
		coe=Roots[iy]
		A=MPSM(ix,coe)
		pushfirst!(MPS, A)
		####################
		if length(MPS)>1
			MPOx=[MPO0[1:iy]; MPO]
		 	MPS,dimc=evolution(MPOx,MPS,sl)
		end
		sl=1-sl
	end
	T=MPOM(ix+1)
	for iy=1:Ol
		append!(MPO,[T])
	end
	######## export the maximum bond dimension
	write(fileE, join([(ix+1)*Dt, dimc], " "), "\n")
	##########################################
	#calculate the propgator and rhot
	if div(ix+1,10)*10==ix+1
		Gt=propagator(MPS)	
		write(file, join(Gt, " "), "\n")
		namet="tmp/GT-"*string(ix+1)*"-dim="*string(dim)*"-alpha="*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*".jld"
		jldopen(namet, "w") do file
			write(file, "ix", ix+1)
			write(file, "dimc", dimc)
			write(file, "Gt", Gt)
		end
		t=(ix+1)*Dt
		###return the Schordinger picture
		rhot=exp(t*Hsn)*Gt*rho0
		#################################
		Sz=[0,0,0,1]'*rhot
		S0=[1,0,0,0]'*rhot
		println("ix=", ix+1, " dimc=", dimc, " Sz=", Sz, " S0=",S0)
	end
end
###############################################
M=Matrix{ComplexF64}(I,4,4)
for ix=1:Nt
	global M
	Tx=MPS[ix]
	M,FD=bondDimension(M,Tx,"l")
	write(fileD,join(FD, " "), "\n")
	#println("BD=",join(FD, " "))
end
###############################################
#stroing the data
#jldopen("data/GTensor-alpha-"*string(alphac)*"-K="*string(K)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-ep="*string(epsilon)*"-HD.jld2","w") do file
#	for i=1:Nt
#		write(file, "G"*string(i),MPS[i])
#	end
#end
##############################################
#Gt=propagator(MPS)
#rhot=reshape(Gt*rho0,(2,2))
#println("ix=", ix, " mxi rhot = ", tr(rhot*rhot))
####exact solution
#	A1=MPSM(0,1.0+0.0im)
#	A2=MPSM(1,1.0+0.0im)
#	tx=0*Dtx
#	Cpp0=alphac*(log(1+2*(tx^2+Dtx^2)+(tx^2-Dtx^2)^2)-2*log(1+tx^2))/(2*Dtx)^2
#	Cpn0=alphac*(2*atan(tx)-atan(tx-Dtx)-atan(tx+Dtx))/Dtx^2
#	GtE=A1[:,1,:]+Cpp0/2.0.*A1[:,2,:]*A1[:,2,:]
#	rhotE=reshape(GtE*rho0,(2,2))
#	println("exact mxi rhot0 = ", tr(rhotE*rhotE))
#	tx=Dtx
#	Cpp1=alphac*(log(1+2*(tx^2+Dtx^2)+(tx^2-Dtx^2)^2)-2*log(1+tx^2))/(2*Dtx)^2
#	Cpn1=alphac*(2*atan(tx)-atan(tx-Dtx)-atan(tx+Dtx))/Dtx^2
#	xpp=1.0/2.0*Cpp1*Cpp1.*A2[:,2,:]*A2[:,2,:]*A1[:,2,:]*A1[:,2,:]
#	xnp=1.0/2.0*Cpp1*Cpn1.*A2[:,2,:]*A2[:,2,:]*A1[:,2,:]*A1[:,3,:]
#	xpn=1.0/2.0*Cpp1*Cpn1.*A2[:,2,:]*A2[:,2,:]*A1[:,3,:]*A1[:,2,:]
#	xnn=1.0/2.0*Cpn1*Cpn1.*A2[:,2,:]*A2[:,2,:]*A1[:,3,:]*A1[:,3,:]
#	GtE=(A2[:,1,:]+Cpp0/2.0.*A2[:,2,:]*A2[:,2,:])*(A1[:,1,:]+Cpp0/2.0.*A1[:,2,:]*A1[:,2,:])+A2[:,2,:]*(Cpp1.*A1[:,2,:]+Cpn1.*A1[:,3,:])+xpp+xpn+xnp+xnn
#	rhotE=reshape(GtE*rho0,(2,2))
#	println("exact mxi rhot1 = ", tr(rhotE*rhotE))
