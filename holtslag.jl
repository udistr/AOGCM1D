include("constants.jl")
include("VelocityScale.jl")
include("holtslag_psim.jl")
include("holtslag_psis.jl")
include("holtslag_local.jl")

function holtslag(θv,U,V,LH,SH,ustar,ρA,ZF)
# a vertical diffusion scheme based on Holtslag and Boville (1993)
# Includes enhenced diffusion below the PBL
#θv-virtual potential temperature

Rcrit=0.5;
b=8.5
eps=0.1;
minshear=0.003
d=10

Z=ZF[2:end];
DZ=diff(ZF,dims=1);
ZM=Z-DZ./2;
SZM=length(ZM);

ZD=[Z[1]:0.01*(Z[1]):Z[end];]

#wθv-surface virtual heat flux, positive upward
wθv=SH./ρA./cpa.+humid_fac.*θv[1].*LH./ρA./Av;

#TF=(SH+LH)./ρA./cpa;
ustar=abs(ustar)
L=-θv[1].*ustar.^3 ./(karman.*gravity_mks.*wθv);
HD=ZD./L;
H=ZF[2:end]./L;

# function to calculate velocity VelocityScale

# surface velocity scale
out=VelocityScale.(0.1*ZD,ZD,L,ustar,θv[1],wθv)
ws=map(x->x[1], out)
wm=map(x->x[2], out)

θs=θv[1].+max.(b.*wθv./ws,0);

# interpolation to high resolution
spl = Spline1D(ZM, U); UD=spl(ZD);
spl = Spline1D(ZM, V);  VD=spl(ZD);
spl = Spline1D(ZM, θv); TD=spl(ZD);
#calculate PBLH at high resolution
ha=(Rcrit.*θs.*(UD.^2+VD.^2)./(gravity_mks.*(TD-θs)))-ZD;
ha[ha.>0].=Inf;
if all(isinf.(ha))
    if wθv>0
        id=SZM;
    else
        id=1;
    end
else
    ~,id=findmin(abs.(ha));
end
id=max(id,1);
h=abs(ZD[id[1]]);
# n - levels inside PBL
n=(h.-ZF).*((h.-ZF).>0);
n[n.==0].=Inf;
~,i=findmin(n);
# frac - fraction of the upper PBL level that is indide the PBL
frac=(h-ZF[i])./DZ[i];

# nsl - levels inside surface layer
nsl=(0.1 .*h.-ZF).*((0.1 .*h.-ZF).>0);
nsl[nsl.==0].=Inf;
~,isl=findmin(nsl);

#ipbl - levels inside the PBL
#npbl - levels outside the PBL
ipbl=[1:i;];
if i==SZM
    npbl=[];
else
    npbl=i+1:SZM;
end

kM=zeros(SZM); kS=zeros(SZM);
wm=zeros(SZM); wt=zeros(SZM);

out=VelocityScale.(ZM[ipbl],h,L,ustar,θv[1],wθv)

wt[ipbl]=map(x->x[1], out)
wm[ipbl]=map(x->x[2], out)
wstar=((gravity_mks/θv[1])*wθv*(wθv>0)*h)^(1/3)

kS[ipbl]=karman.*wt[ipbl].*min.(ZM[ipbl],h).*(1 .- min.(ZM[ipbl]/h,1)).^2
kM[ipbl]=karman.*wm[ipbl].*min.(ZM[ipbl],h).*(1 .- min.(ZM[ipbl]/h,1)).^2

kM[i]=kM[i].*frac
kS[i]=kS[i].*frac

#vertical diffusion above the PBL
S=max.(abs.(diff([0; U])./diff(ZF)),minshear);
GRi=(gravity_mks./θv[1:end-1]).*diff(θv)./diff(ZM)./S[1:end-1].^2;

fm=[0; holtslag_local.(GRi)]

lc=30

i2=min(i+1,SZM)
kM[npbl]=kM[npbl]+lc.^2 .*S[npbl].*fm[npbl];
kS[npbl]=kM[npbl]+lc.^2 .*S[npbl].*fm[npbl];

a=7.2

γcq=a*(LH./ρA./Av)/wm.^2*wstar/h
γct=a*(SH./ρA./cpa)/wm.^2*wstar/h
γcm=a*(ustar^2)/wstar*wm.^2/h


println(h)
return kM,kS,γcq,γct,γcm
#[i max(kM) max(kS) wθv max(Ts) U(1)]

end

  