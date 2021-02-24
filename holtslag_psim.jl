function holtslag_psim(xi)
if (xi>=0) & (xi<=1.)
    y=1+5*xi;
elseif xi>=1
    y=5+xi
else xi<0
    y=(1-15*xi)^(-1/3)
end
return y
end
