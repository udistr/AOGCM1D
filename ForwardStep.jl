using Dierckx
using Plots
using RollingFunctions
using LinearAlgebra
using Dierckx
using Plots
using RollingFunctions
using LinearAlgebra
using Dates

function ForwardStep(mtime)

  global Θv,TA,SST,ρA,ρO,ustar,LH,SH,E,LW,TAU,Qnet,SW
  global KAm,KAt,KOm,KOt
  global UA1,qA1,ΘA1,UO1,TO1
  global qA,ΘA,UA,UO,TO,SO

  ########################################################################
  # time
  ########################################################################

  epochdays=Dates.date2epochdays(Date(mtime))
  t=(Dates.datetime2epochms(mtime)/1000-epochdays*24*60*60)/3600;
  #julian day
  jd=Dates.datetime2julian(mtime)-Dates.datetime2julian(DateTime(year(mtime)));

  ########################################################################
  # radiation
  ########################################################################
  AOI=angle_of_incidence(lat,lon,jd,t); 
  AOI<0 ? AOI=0 : AOI
  f=orad(ZOF)
  SW=(SunConstant*(1-albedo).*f.*AOI);

  ########################################################################
  # surface
  ########################################################################

  hl,hs,_,_,_,_,_,_,_,_,_,ustar,_,_,_,_=
  bulk(TA[1]+d2k,qA[1],abs(UA[1]-UO[1]),
        SST,ZAC[1],ZAC[1],ZAC[1],ZAC[1],ρA[1]);
  LH=-hl;
  SH=-hs;
  E=-hl/Av;
  ϵa=0.725+0.17*log10(sum(qA/1000*1.2*ΔZA)*100)
  LW=emissivity*sb*((SST+d2k).^4 .- ϵa.*(TAup).^4);
  TAU=ρA[1].*ustar.^2 .*sign(ustar);
  Qnet=SW[1]-LH-SH-LW;
  println(Qnet)

  ########################################################################
  # Atmospheric vertical diffusion
  ########################################################################
  Θv=ΘA.*(1 .+humid_fac*qA);
  KAm,KAt,γcq,γct,γcm=holtslag(Θv,UA,UA.*0,LH,SH,ustar,ρA[1],ZAF);
  KAm,KAt,γcq,γct,γcm=holtslag(Θv,UA,UA.*0,LH,SH,ustar,ρA[1],ZAF);
  #println("----",maximum(γct),"----")
  ########################################################################
  # Oceanic vertical diffusion
  ########################################################################
  ρO=ρO0.*(1 .-alpha.*(TO.-tref).+beta.*(SO.-sref));
  B=gravity_mks.*(ρO0.-ρO)./ρO0;
  KOm,KOt,_= large(TO,TO.*0,UO,VO,ρO,Qnet,0,ustar,ZOF,lat);

  ########################################################################
  # Atmospheric explicit time stepping
  ########################################################################
  TAUS=sign(UA[1]-UO[1]);

  UAst[1]=UA[1]-TAU*TAUS/(ΔZA*ρA[1])*ΔT-(1-Aimp)*(
      KAm[1].*(UA[1]-UA[2])/ΔZA^2*ΔT);
  UAst[end]=UA[end]-(1-Aimp)*(
      KAm[end]*(UA[end]-UAup)/ΔZA^2*ΔT);
  UAst[2:end-1]=UA[2:end-1]+(1-Aimp)*(
      KAm[1:end-2].*(UA[1:end-2]-UA[2:end-1])/ΔZA^2*ΔT-
      KAm[2:end-1].*(UA[2:end-1]-UA[3:end  ])/ΔZA^2*ΔT);

  qAst[1]=qA[1]-E/(ΔZA*ρO[1])*ΔT-
      (1-Aimp)*(KAt[1]*(qA[1]-qA[2])/ΔZA^2*ΔT+
      (-γcq[1].*KAt[1])/ΔZA*ΔT);
  qAst[end]=qA[end]-(1-Aimp)*(KAt[end]*(qA[end]-q0)/ΔZA^2*ΔT+
      (γcq[end-1].*KAt[end-1]-γcq[end].*KAt[end])/ΔZA*ΔT);
  qAst[2:end-1]=qA[2:end-1]+(1-Aimp)*(
      KAt[1:end-2].*(qA[1:end-2]-qA[2:end-1])/ΔZA^2*ΔT-
      KAt[2:end-1].*(qA[2:end-1]-qA[3:end  ])/ΔZA^2*ΔT+
      (γcq[1:end-2].*KAt[1:end-2]-γcq[2:end-1].*KAt[2:end-1])/ΔZA*ΔT);

  ΘAst[1]=ΘA[1]+(SH)/cpa/(ΔZA*qA[1])*ΔT-
        (1-Aimp)*(KAt[1]*(ΘA[1]-ΘA[2])/ΔZA^2*ΔT+
        (-γct[1].*KAt[1])/ΔZA*ΔT);
  ΘAst[end]=ΘA[end]-
        (1-Aimp)*(KAt[end]*(ΘA[end]-ΘA0)/ΔZA^2*ΔT+
        (γct[end-1].*KAt[end-1]-γct[end].*KAt[end])/ΔZA*ΔT);
  ΘAst[2:end-1]=ΘA[2:end-1]+(1-Aimp)*(
        KAt[1:end-2].*(ΘA[1:end-2]-ΘA[2:end-1])/ΔZA^2*ΔT-
        KAt[2:end-1].*(ΘA[2:end-1]-ΘA[3:end  ])/ΔZA^2*ΔT+
        (γct[1:end-2].*KAt[1:end-2]-γct[2:end-1].*KAt[2:end-1])/ΔZA*ΔT);

  ########################################################################
  # Oceanic explicit time stepping
  ########################################################################
  if cpl==1

      TOst[1]=TO[1]+(SW[1]-LH-SH-LW)/cp/(ΔZO*ρO[1])*ΔT-
          (1-Oimp)*(KOt[1]*(TO[1]-TO[2])/ΔZO^2*ΔT);
      TOst[end]=TO[end]+SW[end]/cp/(ΔZO*ρO[end])*ΔT+
          (1-Oimp)*(KOt[end-1].*(TO[end-1]-TO[end])/ΔZO^2*ΔT);
      TOst[2:end-1]=TO[2:end-1]+(1-Oimp)*(
          KOt[1:end-2].*(TO[1:end-2]-TO[2:end-1])/ΔZO^2*ΔT-
          KOt[2:end-1].*(TO[2:end-1]-TO[3:end  ])/ΔZO^2*ΔT)+
          SW[2:end-1]./cp./(ΔZO*ρO[2:end-1])*ΔT;

      UOst[1]=UO[1]+
            (TAU)*TAUS/(ρO[1])/ΔZO*ΔT-
            (1-Oimp)*(KOm[1]*(UO[1]-UO[2])/ΔZO^2*ΔT);
      UOst[end]=UO[end]+
            (1-Oimp)*(KOm[end-1].*(UO[end-1]-UO[end])/ΔZO^2*ΔT);
      UOst[2:end-1]=UO[2:end-1]+(1-Oimp)*(
            KOm[1:end-2].*(UO[1:end-2]-UO[2:end-1])/ΔZO^2*ΔT-
            KOm[2:end-1].*(UO[2:end-1]-UO[3:end  ])/ΔZO^2*ΔT);

  end

  if Aimp==0
      UA1=UAst; qA1=qAst; ΘA1=ΘAst;
  else
      du=-Aimp*KAm[1:end-1].*ΔT./ΔZA^2; dl=-Aimp*KAm[2:end].*ΔT./ΔZA^2;
      d=(1 .+Aimp.*ΔT./ΔZA^2 .*([KAm[1:end-1]; 0]+[0; KAm[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      UA1=A*UAst;
      du=-Aimp*KAt[1:end-1].*ΔT./ΔZA^2; dl=-Aimp*KAt[2:end].*ΔT./ΔZA^2;
      d=(1 .+Aimp.*ΔT./ΔZA^2 .*([KAt[1:end-1]; 0]+[0; KAt[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      qA1=A*qAst; ΘA1=A*ΘAst;
  end

  if Oimp==0
      UO1=UOst; TO1=TOst;
  else
      du=-Oimp*KOm[1:end-1].*ΔT./ΔZO^2; dl=-Oimp*KOm[2:end].*ΔT./ΔZO^2;
      d=(1 .+Oimp.*ΔT./ΔZO^2 .*([KOm[1:end-1]; 0]+[0; KOm[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      UO1=A*UOst;
      du=-Oimp*KOt[1:end-1].*ΔT./ΔZO^2; dl=-Oimp*KOt[2:end].*ΔT./ΔZO^2;
      d=(1 .+Oimp.*ΔT./ΔZO^2 .*([KOt[1:end-1]; 0]+[0; KOt[2:end]]));
      A=inv(Tridiagonal(dl,d,du));
      TO1=A*TOst;
  end
  SST=TO[1]

  TA=(ΘA1).*(PA./PA[1]).^(2/7).-d2k;
  ρA=PA./(Rgas.*(TA.+d2k));

  if moist==1
      qsat=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA.+d2k))./ρA;
      dq=(qA.-qsat).*((qA1.-qsat).>0);
      dq=min.(dq,qA1);
      dTHETA=Av.*dq./cp; #[J/kg][j-1 kg K]
      #println([i maximum(dTHETA) maximum(dq)])
      ind=qA1.>qsat;
      qA1[ind]=qA1[ind]-dq[ind];
      #THETA[i,:]=THETA[i,: ]+dTHETA;
  end

  qA=qA1; ΘA=ΘA1; UA=UA1;
  UO=UO1; TO=TO1;
end