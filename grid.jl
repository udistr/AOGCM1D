# grid
lat=30
lon=0

# atmospheric grid 
ΔZA=50; # atmospheric layer thickness
vleva=50; # number of atmospheric vertical levels
ZAF=ΔZA*[0:(vleva);];
ZAC=ZAF[1:end-1].+ΔZA/2

# ocean grid 
ΔZO=1; # oceanic layer thickness
vlevo=100; # number of oceanic vertical levels
# edges
ZOF=zeros(vlevo+1);
ZOF[2:vlevo+1]=-[1:vlevo;]*ΔZO;
ZOC=ZOF[1:end-1].+ΔZO/2