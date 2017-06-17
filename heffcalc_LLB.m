%calculation of effective field and coefficients in LLB equation
%input
%1. D,[J] double, uniaxial anisotropy
%2. muRE, muTM, [J/T], double, magnetic moment for RE(rare earth), TM(transition metal)
%3. mRE, mTM, [unitless], 1x3 vector, magnetization of RE and TM
%4. mTTM, mTRE,[unitless], 1x3 vector, equilibrium magnetization of RE and
%TM, calculated from curie weiss equation
%5. Hext, [Tesla], 1x3 vector, external field
%8. MsT, [A/m],double, total magnetic moment at zero temperature
%9. ita,[unitless],double, torque efficiency
%10. PFL,[unitless],double, FL polarization
%11. Jc,[A/m2],double, charge current density pass through HM
%12. elev,[eV],double, electron charge
%13. tFL,[m],double, FL thickness
%14. J0RERE,J0TMTM,J0RETM,J0TMRE, [J],double, exchange coefficient
%15. alp,[unitless],double, damping constant
%16. ip,[unitless],1x3 vector, spin current polarization
%17. bbeta,[J],double,1/kbT
%18. addSTT, [1/0], flag if include STT
%19. addSOT, same as addSTT
%20. thetaSHE, [unitless], double, spin hall angle, 
%output
function [Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA0,...
    HRE_MFA0,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,MsT,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta,addSTT,addSOT,thetaSHE)
constantfile();
HARE=[0,0,2*D/muRE*mRE(3)];
HATM=[0,0,2*D/muTM*mTM(3)];
HeffRE=Hext+HARE;
HeffTM=Hext+HATM;

if size(MsT,2)==1%for scalar (mz) solution of curie weiss equation
   if addSOT
    Js=thetaSHE*Jc;
Hi=Js*hbar/(2*elev*MsT*tFL);   
elseif addSTT
Hi=ita*PFL*Jc*hbar/(2*elev*MsT*tFL);
   end
else %for vector solution of curie weiss equation
    if mTRE<1e-3 && mTTM<1e-3
    error('Ms=0,torque singularity')
else
%MsT=abs((x*muRE/mub*mTRE+q*muTM/mub*mTTM))/Msperatom*Ms0;
end
Hitmp=zeros(1,3);
for cttmp=1:size(mTTM,2)
    if MsT(cttmp)==0
        Hitmp(cttmp)=0;
    else
        Hitmp(cttmp)=ita*PFL*Jc*hbar./(2*elev*MsT(cttmp)*tFL);
    end
end
Hi=[Hitmp(1),Hitmp(2),Hitmp(3)];
end
HRE_MFA=(muRE*HeffRE+J0RERE*mRE+J0RETM*mTM+muRE*Hi.*ip/alp)/muRE;%eqn5
HTM_MFA=(muTM*HeffTM+J0TMTM*mTM+J0TMRE*mRE+muTM*Hi.*ip/alp)/muTM;%eqn6
HRE_MFA0=(muRE*HeffRE+J0RERE*mRE+J0RETM*mTM)/muRE;%eqn5
HTM_MFA0=(muTM*HeffTM+J0TMTM*mTM+J0TMRE*mRE)/muTM;%eqn6

xi0_TM=bbeta*muTM*HTM_MFA;%eqn5
xi0_RE=bbeta*muRE*HRE_MFA;
xi0_TM0=bbeta*muTM*HTM_MFA0;%eqn5
xi0_RE0=bbeta*muRE*HRE_MFA0;
% [Bri_val_TM,Bri_prime_TM]=Bri_func(norm(xi0_TM),JFe);
% [Bri_val_RE,Bri_prime_RE]=Bri_func(norm(xi0_RE),JGd);
[Bri_val_TM,~]=Lang_fun(norm(xi0_TM));
[Bri_val_RE,~]=Lang_fun(norm(xi0_RE));
[Bri_val_TM0,Bri_prime_TM0]=Lang_fun(norm(xi0_TM0));
[Bri_val_RE0,Bri_prime_RE0]=Lang_fun(norm(xi0_RE0));

m0_TM=Bri_val_TM*xi0_TM/norm(xi0_TM);%eqn5
m0_RE=Bri_val_RE*xi0_RE/norm(xi0_RE);
% m0_TM0=Bri_val_TM0*xi0_TM0/norm(xi0_TM0);%eqn5
% m0_RE0=Bri_val_RE0*xi0_RE0/norm(xi0_RE0);

Lambda_TM=2*gam*alp/(bbeta*muTM);%eqn.6
Lambda_RE=2*gam*alp/(bbeta*muRE);

Gam_parall_TM=Lambda_TM*Bri_val_TM0/(norm(xi0_TM0)*Bri_prime_TM0);
Gam_parall_RE=Lambda_RE*Bri_val_RE0/(norm(xi0_RE0)*Bri_prime_RE0);

Gam_perp_TM=Lambda_TM/2*(norm(xi0_TM0)/Bri_val_TM0-1);
Gam_perp_RE=Lambda_RE/2*(norm(xi0_RE0)/Bri_val_RE0-1);
end