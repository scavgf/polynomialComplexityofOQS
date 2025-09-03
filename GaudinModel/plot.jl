include("./parameter.jl")
using .Parameters:Omega, Dt,Nt,dim,N,A,epsilon,Od,lambda
using Plots
using LaTeXStrings


k=30
fileBD=open("data/propagator-Omega="*string(2*Omega)*"-N="*string(k)*"-D-"*string(dim)*"-Nt="*string(Nt)*"-Dt="*string(Dt)*"-lambda="*string(lambda)*"-ep="*string(epsilon)*"-Od="*string(Od)*".dat","r")
data=readlines(fileBD)
close(fileBD)
xlist=Vector{Float64}()
ylist=Vector{Float64}()
for ix=1:length(data)   
    Gt=[parse(Float64,x) for x in split(data[ix], " ")]
    Szt=last(Gt)/2.0
    append!(xlist, ix*Dt)
    append!(ylist, Szt)
end
plot(xlist, ylist)
xlabel!(L"J t")
ylabel!(L"S_z")
###close the file



