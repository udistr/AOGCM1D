include("holtslag_psim.jl")
include("holtslag_psis.jl")
include("constants.jl")

function VelocityScale(z,h,L,ustar,θv=nothing,wθv=nothing)
  xi=z/L
  psis=holtslag_psis(xi)
  psim=holtslag_psim(xi)
  if xi>0
      wt=ustar/psis
      wm=wt
  else
      if z/h<0.1
          wt=ustar/psis
          wm=ustar/psim
      else
          c1=0.6
          a=7.2
          wθv<0 ? println("wθv",wθv) : nothing
          #println("θv",θv)
          wstar=((gravity_mks/θv)*wθv*h)^(1/3)
          wm=(ustar^3+c1*wstar^3)^(1/3)
          xi1=0.1*z/L
          psis1=holtslag_psis(xi1)
          psim1=holtslag_psim(xi1)
          Pr=psis1/psim1*abs(xi1)+a*karman*(0.1*(wstar/wm))
          wt=wm/Pr
      end
  end    
  return wt,wm
end