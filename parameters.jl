include("constants.jl")

# schemes
cpl=1;
moist=1;

#boundary conditions
TAup=270; #top model temperature
UAup=(0.2/karman).*(log((ZAF[end]+ΔZA[end]/2)/1e-4)); #top model temperature
UOdown=0; #lower boundary model current

# numerics
Aimp=0.5; # Atmospheric implicitness
Oimp=0.5; # Oceanic implicitness