using Dierckx
using Plots
using RollingFunctions
using LinearAlgebra
using Dates

include("constants.jl")
include("grid.jl")
include("parameters.jl")
include("time.jl")

include("delta.jl")
include("heaviside.jl")
include("angle_of_incidence.jl")
include("bulk.jl")
include("holtslag.jl")
include("large.jl")
include("orad.jl")
include("initialize.jl")
include("ForwardStep.jl")

mtime=stime # model clock
Q=zeros(n)
sw=zeros(n)
lw=zeros(n)
lh=zeros(n)
sh=zeros(n)
for i=2:n+1

    global mtime
    println(i,"   ",mtime)
    ForwardStep(mtime)
    
    if rem(i,360)==0;#3600*6/dt
      title = plot(title = mtime, grid = false, axis = false,
       showaxis = false, bottom_margin = -50Plots.px,yaxis=nothing)
      f1=plot(-ZOC,TO,ylabel=("TO"));#,ylims=(19,21)
      f2=plot(-ZOC,UO,ylabel=("UO"));
      f3=plot(-ZOC,KOm,ylabel=("KOm"));
      f4=plot(ZAC,UA,ylabel=("UA"));
      f5=plot(ZAC,Θv.-273.15,ylabel=("Θv")); #ylim([19 21]);...
      f6=plot(ZAC,qA,ylabel=("q"));
      f7=plot(ZAC,KAm,ylabel=("KAm"));
      plt=plot(title,f1,f2,f3,f4,f5,f6,f7,layout=(8,1),legend=false,titlefontsize=6,
        size = (500, 1000))
      display(plt)
    end
    mtime=mtime+Second(ΔT)
    Q[i-1]=Qnet
    sw[i-1]=sum(SW)
    lw[i-1]=LW
    lh[i-1]=LH
    sh[i-1]=SH
end