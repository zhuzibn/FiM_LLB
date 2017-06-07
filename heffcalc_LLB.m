%calculation of effective field and coefficients in LLB equation
%input
%1. D,[J] double, uniaxial anisotropy
%2. muRE, muTM, [J/T], double, magnetic moment for RE(rare earth), TM(transition metal)
%3. mRE, mTM, [unitless], 1x3 vector, magnetization of RE and TM
%4. mTTM, mTRE,[unitless], 1x3 vector, equilibrium magnetization of RE and
%TM, calculated from curie weiss equation
%5. Hext, [Tesla], 1x3 vector, external field 
%6. x,q, [unitless], double, concentration of RE,TM
%7. Msperatom,[mub],double, total magnetic moment at zero temperature
%8. Ms0, [A/m],double, total magnetic moment at zero temperature
%9. ita,[unitless],double, torque efficiency
%10. PFL,[unitless],double, FL polarization
%11. Jc,[A/m2],double, charge current density pass through HM
%12. elev,[eV],double, electron charge
%13. tFL,[m],double, FL thickness
%14. J0RERE,J0TMTM,J0RETM,J0TMRE, [J],double, exchange coefficient
%15. alp,[unitless],double, damping constant
%16. ip,[unitless],1x3 vector, spin current polarization
%17. bbeta,[J],double,1/kbT
%output
function [Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA,...
    HRE_MFA,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,x,q,Msperatom,Ms0,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta)
constantfile();
kb=1.38064852e-23;%[J/K]
mub=9.274e-24;%[J/T]bohr magneton
hbar=6.58211951440e-16;%[eV.s]
elev=1;%[electron charge]
gam=1.760859644e11;%[rad/(s.T)]
HARE=[0,0,2*D/muRE*mRE(3)];
HATM=[0,0,2*D/muTM*mTM(3)];
HeffRE=Hext+HARE;
HeffTM=Hext+HATM;

MsT=abs((x*muRE/mub*mTRE+q*muTM/mub*mTTM))/Msperatom*Ms0;  
Hi=ita*PFL*Jc*hbar/(2*elev*MsT*tFL);

HRE_MFA=(muRE*HeffRE+J0RERE*mRE+J0RETM*mTM+muRE*Hi/alp*ip)/muRE;%eqn5
HTM_MFA=(muTM*HeffTM+J0TMTM*mTM+J0TMRE*mRE+muTM*Hi/alp*ip)/muTM;%eqn6

xi0_TM=bbeta*muTM*HTM_MFA;%eqn5
xi0_RE=bbeta*muRE*HRE_MFA;

% [Bri_val_TM,Bri_prime_TM]=Bri_func(norm(xi0_TM),JFe);
% [Bri_val_RE,Bri_prime_RE]=Bri_func(norm(xi0_RE),JGd);
[Bri_val_TM,Bri_prime_TM]=Lang_fun(norm(xi0_TM));
[Bri_val_RE,Bri_prime_RE]=Lang_fun(norm(xi0_RE));


m0_TM=Bri_val_TM*xi0_TM/norm(xi0_TM);%eqn5
m0_RE=Bri_val_RE*xi0_RE/norm(xi0_RE);

Lambda_TM=2*gam*alp/(bbeta*muTM);%eqn.6
Lambda_RE=2*gam*alp/(bbeta*muRE);

Gam_parall_TM=Lambda_TM*Bri_val_TM/(norm(xi0_TM)*Bri_prime_TM);
Gam_parall_RE=Lambda_RE*Bri_val_RE/(norm(xi0_RE)*Bri_prime_RE);

Gam_perp_TM=Lambda_TM/2*(norm(xi0_TM)/Bri_val_TM-1);
Gam_perp_RE=Lambda_RE/2*(norm(xi0_RE)/Bri_val_RE-1);
end