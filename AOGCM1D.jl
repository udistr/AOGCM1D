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
#=
k1=plot(rollmean(U0[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="U0");
k2=plot(rollmean(UA[:,end],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="Utop");
k3=plot(rollmean(UA[:,end]-U0[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="Utop-U0");
k4=plot(rollmean(UA[:,1],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="U");
k5=plot(rollmean(THETA[:,2],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="THETA");
k6=plot(rollmean(SST,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="SST");
k7=plot(rollmean(TO[:,10],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="T");
#tt=string("ΔZO=",ΔZO)
plt=plot(k1,k2,k3,k4,k5,k6,k7,layout=(7,1),legend=false)
display(plt)

#PyPlot.suptitle(tt)

k1=plot(rollmean(xi[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="xi");
k2=plot(rollmean(ustar[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="ustar");
k3=plot(rollmean(psimh[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="psimh");
k4=plot(rollmean(psixh[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="psixh");
k5=plot(rollmean((rd[1,:].^2),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CD");
k6=plot(rollmean((rh[1,:].*rd[1,:]),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CH");
k7=plot(rollmean((re[1,:].*rd[1,:]),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CE");
#tt=string("ΔZO=",ΔZO)
plt=plot(k1,k2,k3,k4,k5,k6,k7,layout=(7,1),legend=false)
display(plt)

#PyPlot.suptitle(tt)

k1=plot(rollmean(SW[:,1]/cp/(ΔZO*ρO0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="SW");
k2=plot(rollmean(LW[1,:]/cp/(ΔZO*ρO0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="LW");
k3=plot(rollmean(SH[1,:]/cp/(ΔZO*ρO0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="SH");
k4=plot(rollmean(LH[1,:]/cp/(ΔZO*ρO0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="LH");
k5=plot(rollmean((SW[:,1]-LW[1,:]-SH[1,:]-LH[1,:])/cp/(ΔZO*ρO0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="QNET");
plt=plot(k1,k2,k3,k4,k5,layout=(5,1),legend=false);
display(plt)

#PyPlot.suptitle(tt)
=#