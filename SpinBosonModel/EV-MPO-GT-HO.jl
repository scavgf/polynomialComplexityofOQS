module Evolution

include("./parameter.jl")
using LinearAlgebra
using TensorOperations 
using TimerOutputs
using .Parameters:dim, epsilon, Nt, dim, Dt, K,alphac

export evolution 
##########################################################################################################
function evolution(MPO::Vector{Any}, MPS::Vector{Any},  dir::Int64)
	#the size of Mr and Mr
	######get the running status
	dimc=0
	L=length(MPO)
	if dir==0
		SL=MPS[1] #the left tensor
		ML=complex(MPO[1]) #the boundary M matrix
		#RG from left to right
		for jx=2:L-1
			M=complex(MPO[jx])
			MU,MV=decomposition(M,"l")
			SR=MPS[jx] #the left tensor 
			MR=MU # the left M matrix
			#RG SL and SR according to MPOs: ML and MR
			SL0,SL,dimx=LRRG(SL,SR,ML,MR,"l")
			#update the right one
			MPS[jx-1]=SL0
			#update MR matrix
			ML=MV
			###output the largest bond-dimension
			dimc= dimc> dimx ? dimc : dimx
		end
		#update the right boundary
		#######################################
		SR=MPS[L]
		MR=complex(MPO[L][:,:,:,1])
		#######################################
		SL0,SL,dimx=LRRG(SL,SR,ML,MR,"l")
		MPS[L-1]=SL0
		MPS[L]=SL
		###output the largest bond-dimension
		dimc= dimc> dimx ? dimc : dimx
		#####################################

	else
		#######################################
		#RG from right to left (Do the complex conjucated propocation)
		#starting from right boundary
		MR=complex(MPO[L][:,:,:,1])
		SR=MPS[L]
		for jx=L-1:-1:2 #update MPS from right to left
			M=complex(MPO[jx])
			MU,MV=decomposition(M,"r")
			ML=MV
			SL=complex(MPS[jx])
			#######################################
			SR,SR0,dimx=LRRG(SL,SR,ML,MR,"r")
			#######################################
			MPS[jx+1]=SR0
			MR=MU
			####determine the lagest bond dimension####
			dimc= dimc> dimx ? dimc : dimx
		end
		#update the left boundary
		SL=MPS[1]
		ML=complex(MPO[1])
		SR,SR0,dimx = LRRG(SL,SR,ML,MR,"r")
		####determine the lagest bond dimension####
		dimc= dimc> dimx ? dimc : dimx
		MPS[2]=SR0
		MPS[1]=SR
	end
	return MPS, dimc
end
##############################################################################################################
function decomposition(M::Array{ComplexF64,4},sgn::String)
	##do a SVD of M from###
	sm=size(M)
	if sgn=="r"
		#from right to left
		Mmr=reshape(M,(sm[1]*sm[2],sm[3]*sm[4]))
		Fr=svd(Mmr)
		Dr=Diagonal(sqrt.(Fr.S))
		MU=reshape(Fr.U*Dr, (sm[1], sm[2], sm[1]*sm[2]))
		MV=reshape(Dr*Fr.Vt, (sm[1]*sm[2], sm[3], sm[4]))
	else
		#from left to right
		Mml=reshape(permutedims(M,(1,3,2,4)),(sm[1]*sm[3],sm[2]*sm[4]))
		Fl=svd(Mml)
		Dl=Diagonal(sqrt.(Fl.S))
		MU=permutedims(reshape(Fl.U*Dl,(sm[1],sm[3],sm[1]*sm[3])),(1,3,2))
		MV=permutedims(reshape(Dl*Fl.Vt,(sm[1]*sm[3],sm[2],sm[4])),(2,1,3))
	end
	return MU, MV
end
##############################################################################################################
function LRRG(SL::Array{ComplexF64,3}, SR::Array{ComplexF64,3},  ML::Array{ComplexF64,3}, MR::Array{ComplexF64,3},sgn::String)
	@tensor TL[sout,g,gb,sin]:=ML[g,x,gb]*SL[sout,x,sin]
	@tensor TR[sout,ga,g,sin]:=MR[ga,g,x]*SR[sout,x,sin]
	@tensor SS[sout,gl,gr,sin]:=TL[sout,gl,x,y]*TR[y,x,gr,sin]
    sT=size(SS)
	diml=prod(sT[1:2])
	dimr=prod(sT[3:4])
	TM=reshape(SS,(diml,dimr))
	#svd decomposiiton
	FT=SVD{ComplexF64, ComplexF64, Matrix{ComplexF64}}
    	try
         	FT=svd(TM,alg=LinearAlgebra.DivideAndConquer())
    	catch e
        	FT=svd(TM,alg=LinearAlgebra.QRIteration())
	end
	DS=FT.S
	epsilon1=epsilon*DS[5]
	dimx=length(DS[DS.>epsilon1])
	dimc=(dim<dimx ? dim : dimx)
	#dimc=(dim<length(DS) ? dim : length(DS))
	TU=FT.U[:,1:dimc]
	Tv=FT.Vt[1:dimc,:]
	FS=FT.S[1:dimc]
	Ds=Diagonal(FS)
	if sgn=="r"
		TMl=TU*Ds
		TMr=Tv
	elseif sgn=="l"
		TMl=TU
		TMr=Ds*Tv
	end
	#update the TL
	TL0=reshape(TMl,(sT[1],sT[2],dimc))
	TR0=reshape(TMr,(dimc,sT[3],sT[4]))
	##renormaliztion
	#value= dimc < min(diml,dimr) ? FT.S[dimc+1] : 0
	return TL0, TR0, dimc
end
################################################################################
function propagator(MPS::Vector{Any})
	L=length(MPS)
	LMt=Matrix{Float64}(I,(4,4))
	for i=1:L
		Gxt=(MPS[i])[:,1,:]
		#println("the size of Gt", size(Gxt))
		LMt=LMt*Gxt
	end
	return real(LMt)
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
