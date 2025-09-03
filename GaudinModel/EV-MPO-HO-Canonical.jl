module Evolution

include("./parameter.jl")
using LinearAlgebra
using TensorOperations 
using TimerOutputs
using .Parameters:dim, epsilon, Nt, dim, Dt,Od,lambda

export evolution 
##########################################################################################################
function evolution(MPS::Vector{Array{ComplexF64,3}}, Tk::Array{Float64,4})
	#the size of Mr and Mr
	######get the running status
	dimc=0
	L=length(MPS)
	TL=zeros(ComplexF64,4,4,4)
	TL[:,1,:]=Matrix{ComplexF64}(I,4,4)
	#RG from left to right
	for jx=1:L-1
		TC=MPS[jx] #the left tensor 
		#RG SL and SR according to MPOs: ML and MR
		TL0,TL=LRRG(TL, TC, Tk)
		#update the right one
		MPS[jx]=TL0
	end
	#update the right boundary
	TC=MPS[L]
	Tkr=Tk[:,:,:,1]
	@tensor TL0[l,c,r]:=TL[l,bx,sx]*TC[sx,sc,r]*Tkr[bx,c,sc]
	MPS[L]=TL0
	##################canonicalization###################
	BD=[] #yields the bond dimension
	Ux=Matrix{ComplexF64}(I,4,4)
	for iy=L:-1:2
		T=MPS[iy]
		@tensor Tx[l,c,r]:=T[l, c, x]*Ux[x, r]
		sm=size(Tx)
		M=reshape(Tx,sm[1],sm[2]*sm[3])
		FT=svd(M)
		######bond dimension truncation#####
		DS=FT.S
		epsilon1=epsilon*DS[5]
		dimx=length(DS[DS.>epsilon1])
		dimc=(dim<dimx ? dim : dimx)
		pushfirst!(BD, dimc)
		#######################################
		TU=FT.U[:,1:dimc]
		FS=FT.S[1:dimc]
		Tv=FT.Vt[1:dimc,:]
		###update the MPS and Ux
		MPS[iy]=reshape(Tv, (dimc, sm[2], sm[3]))
		Ux=TU*Diagonal(FS)
	end
	T=MPS[1]
	@tensor Tx[l,c,r]:=T[l,c,x]*Ux[x,r]
	MPS[1]=Tx
	#######################################
	return MPS, BD
end
##############################################################################################################
function LRRG(TL::Array{ComplexF64,3}, TC::Array{ComplexF64,3},  Tk::Array{Float64,4})
	@tensor TT[l,u,r1,r2]:=TL[l,bx,sx]*TC[sx,sc,r2]*Tk[bx,u,sc,r1]
    sT=size(TT)
	diml=prod(sT[1:2])
	dimr=prod(sT[3:4])
	TM=reshape(TT,(diml,dimr))
	#svd decomposiiton
	FT=SVD{ComplexF64, ComplexF64, Matrix{ComplexF64}}
    	try
         	FT=svd(TM,alg=LinearAlgebra.DivideAndConquer())
    	catch e
        	FT=svd(TM,alg=LinearAlgebra.QRIteration())
	end
	DS=FT.S
	epsilon1=epsilon*DS[5]
	dimc=length(DS[DS.>epsilon1])
	#dimc=(dim<length(DS) ? dim : length(DS))
	TU=FT.U[:,1:dimc]
	Tv=FT.Vt[1:dimc,:]
	FS=FT.S[1:dimc]
	Ds=Diagonal(FS)
	TMl=TU
	TMr=Ds*Tv
	#update the TL
	TL0=reshape(TMl,(sT[1],sT[2],dimc))
	TR0=reshape(TMr,(dimc,sT[3],sT[4]))
	##renormaliztion
	#value= dimc < min(diml,dimr) ? FT.S[dimc+1] : 0
	return TL0, TR0
end
################################################################################
function propagator(MPS::Vector{Array{ComplexF64,3}}, coe::ComplexF64)
	##stroe the projector
	Proj=[]
	proj=[2 0 0 0] #trace operation
	for ix=1:Nt-1
		for iy=1:Od-1
			M=MPS[(ix-1)*Od+iy][:,1,:]
			proj=proj*M
		end
			append!(Proj,[proj])
			M=MPS[ix*Od][:,1,:]
			proj=proj*M
	end
	##stroe the propagator
	Prop=[]
	prop=Matrix{ComplexF64}(I,(4,4))
	for ix=Nt:-1:2
		for iy=0:Od-1
			M=MPS[ix*Od-iy][:,1,:]
			prop=M*prop
		end
		###obtain the propgator
		proj=Proj[ix-1]
		M=MPS[(ix-1)*Od][:,:,:]
		Gt=zeros(ComplexF64,4,4)
		Gt[1,:]=proj*M[:,1,:]*prop
		Gt[2,:]=proj*M[:,2,:]*prop./(2*coe*lambda*Dt)
		Gt[3,:]=proj*M[:,3,:]*prop./(2*coe*lambda*Dt)
		Gt[4,:]=proj*M[:,4,:]*prop./(2*coe*lambda*Dt)
		append!(Prop,[Gt])
	end
	return Prop
end
################################################################################
function bondDimension(M::Array{ComplexF64,2},T::Array{ComplexF64,3},sgn::String)
	if sgn=="l"
		@tensor TM[l,c,r]:=M[l,x]*T[x,c,r]
		dim=size(TM)
		MR=reshape(TM,dim[1]*dim[2],dim[3])
	else
		@tensor TM[l,c,r]:=T[l,c,x]*M[x,r]
		dim=size(TM)
		MR=reshape(TM,dim[1],dim[2]*dim[3])
	end
	F=svd(MR)
	if sgn=="l"
		M=F.Vt
	else
		M=F.U
	end
	return M, F.S
end
################################################################################
#end of the module
end
##########################################################################################################
