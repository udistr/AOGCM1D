include("parameters.jl")

############################################################################
# Initial arrays
############################################################################

global TO=zeros(vlevo); global TO1=zeros(vlevo); global TOst=zeros(vlevo);
global UO=zeros(vlevo); global UO1=1zeros(vlevo); global UOst=zeros(vlevo);
global ρO=zeros(vleva);

global UA=zeros(vleva); global UA1=zeros(vleva); global UAst=zeros(vleva);
global qA=zeros(vleva); global qA1=zeros(vleva); global qAst=zeros(vleva);
global ΘA=zeros(vleva); global ΘA1=zeros(vleva); global ΘAst=zeros(vleva);
global ρA=zeros(vleva);
global TA=zeros(vleva);
global PA=zeros(vleva);

############################################################################
# Initial conditions
############################################################################
# atmospheric temperature in DegC
TA=25 .-0.007.*ZAC;
PA[1]=101300;
PA[2:vleva]=PA[1].*exp.(-cumsum(ΔZA./(Rgas.*(TA[2:end].+d2k)./gravity_mks)));
ρA=PA./(Rgas.*(TA.+d2k));
# atmospheric potential temperature in kelvin
ΘA[:]=(TA.+d2k).*(PA[1]./PA).^(2/7);
#ΘA[1:20].=ΘA[20];
ΘA0=ΘA[end].+1;
T0=ones(vlevo).*32;#repmat((16+(vlev-1:-1:0)/vlev*8),n+1,1);
TO=T0[1,1].+15 .*([vlevo:-1:1;]/vlevo.-1)/5;
SST=TO[1];
ssq = saltsat*cvapor_fac*exp(-cvapor_exp./(TA[1]+d2k))./ρA[1];
RH=0.8;
q0=0.5*ssq;
qA=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA.+d2k))./ρA.*RH;
UA=(0.2/karman).*(log.(ZAC./1e-4));
SO=sref

UO=(vlevo:-1:1)/vlevo/10; VO=zeros(size(UO));
atemp=TA[1]+d2k; aqh=qA[1];
speed=abs(UA[1]); sst=SST[1];
tmp0=bulk(atemp,aqh,speed,sst,ZAC[1],ZAC[1],ZAC[1],ZAC[1],ρA[1]);
huol0=tmp0[8];


############################################################################

global LH=zeros(1); global SH=zeros(1); global TAU=zeros(1); TAUO=zeros(1);
global LW=zeros(1); global E=zeros(1); ssq=zeros(1);
xi=zeros(1); ce=zeros(1); ch=zeros(1);
global ustar=zeros(1); qstar=zeros(1); tstar=zeros(1);
psimh=zeros(1); psixh=zeros(1); Ri=zeros(vleva);
rd=zeros(1); rh=zeros(1); re=zeros(1);
global KAm=zeros(vleva); global KAt=zeros(vleva);
global KOm=zeros(vlevo); KOt=zeros(vlevo);
global Qnet=zeros(1);
FWflux=zeros(size(Qnet));