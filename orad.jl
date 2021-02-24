function orad(ZO)

  R=0.62; g1=0.6; g2=20;

  f=R*exp.(ZO/g1)+(1-R)*exp.(ZO/g2);
  f[ZO.<-200].=0;
  f=-diff(f,dims=1);

  return(f)

end
