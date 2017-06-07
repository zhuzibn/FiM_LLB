%sample code for FiM(or FM) with STT and SOT using LLB equation.
%author:ZHU Zhifeng, 2017.June.07
clear all;clc;close all
%[1]. 2012-The Landau-Lifshitz-Bloch equation for ferrimagnetic materials-PRB-U. Atxitia
%[2]. 2009-Magnetic dynamics with spin-transfer torques near the Curie temperature-Paul M. Haney-PRB
%[3]. 2011-Crystallographically amorphous ferrimagnetic alloys Comparing a localized atomistic-PRB-Thomas A. Ostler
%[4]. 1963-Magnetization and electrical resistivity of gadolinium single cry-Harold Eugene Nigh
%[5]. http://mxp.physics.umn.edu/s03/projects/s03ovsm/
%[6]. EE6438, lecture 2
tic
%configuration
%--------------------
mREinitt=1;%1\mRE_init_z=1 0\mRE_init_z=-1
addSTT=0;
addSOT=1;%0\STT model 1\SOT model
if addSTT && addSOT
    error('currently only support STT or SOT');
end
if addSOT
    thetaSHE=0.2;
else
    thetaSHE=0;
end
debbbug=1;
calcFM=0;%calculate property of pure FM
if debbbug+calcFM==0
    bbatch=1;
else
    bbatch=0;
end
if ~(debbbug+calcFM+bbatch==1)
    error('setting error')
end
%--------------------
if debbbug || calcFM
    ppc=2;%\1 7RXR622 2\1PN5LG2 3\Landauer
else
    ppc=3;
end
switch ppc
    case 1
        addpath('D:\dropbox\Dropbox\phd\code\mean_field_approximation');
        addpath('D:\dropbox\Dropbox\phd\code\papers_books\1971-The Magnetization of Pure Iron and Nickel-J. Crangle');
        addpath('D:\dropbox\Dropbox\phd\code\papers_books\1963-Magnetization and electrical resistivity of gadolinium single cry-Harold Eugene Nigh');
    case 2
        addpath('E:\dropbox\Dropbox\phd\code\general\gitcontrol\constant');
        addpath('E:\dropbox\Dropbox\phd\code\general\gitcontrol\mean_field_approximation');
        addpath('E:\dropbox\Dropbox\phd\code\papers_books\1971-The Magnetization of Pure Iron and Nickel-J. Crangle');
        addpath('E:\dropbox\Dropbox\phd\code\papers_books\1963-Magnetization and electrical resistivity of gadolinium single cry-Harold Eugene Nigh');
    case 3
        addpath('/home/a0132576/code/general/mean_field_approximation');
end

if calcFM
    calcTM=1;
    if calcTM
        x_=0;
        T_=[10:40:890,900:10:1200];
    else
        x_=1;
        T_=[10:10:270,272:2:310,312:10:600];
    end
    lang_or_bri=1;
    Hextz_=0;
    runt=500e-12;%[s]
elseif debbbug
    x_=0.25;
    T_=75;%[K]
    lang_or_bri=1;
    if addSTT || addSOT
        Hextz_=0;
        Jc_=50e10;
    else
        Hextz_=-2;
        Jc_=0;
    end
    runt=100e-12;
else
    x_=0.25;
    T_=75;%[K]
    lang_or_bri=1;
    if addSTT
        Hextz_=0;
        Jc_=[1:2:50]*1e10;%[A/m2]
    else
        Jc_=0;
        Hextz_=[-20:1:-8]; %sweep H to get MH loop
    end
    runt=500e-12;%[s]
end

if lang_or_bri
    Jat_TMTM=1.5e-21;Jat_RERE=0.98e-21;
else
    Jat_TMTM=4.5e-21;Jat_RERE=1.26e-21;
end

if calcFM
    Jc_=0;
end
tstep=50e-15;%[s]

szx=size(x_,2);
% constants
constantfile();
kb=1.38064852e-23;%[J/K]
mub=9.274e-24;%[J/T]bohr magneton
hbar=6.58211951440e-16;%[eV.s]
elev=1;%[electron charge]
gam=1.760859644e11;%[rad/(s.T)]
% params
natom=0.001;
z=7;

D=natom*8.07246e-24;%[J]uniaxial anisotropy
JFe=1/2;%p103 of Coey's book
JGd=7/2;%p67 (fig.14) of [4]
JatTMRE=-1.09e-21;
muRE=7.63*mub;muTM=2.217*mub;

alp=0.1;%damping constant, assume same for both element
if addSOT
    %ip=[0,-1,0];
    ip=[0,0,-1];
else
    ip=[0,0,-1];%spin current polarization
end
PFL=0.5;%FL polarization
ita=0.4;%assume fixed torque efficiency
Ms0_Fe=1.71*1e6;%[A/m] from [6]
Ms0_Gd=2639.42*1e3;%[A/m][5]
tFL=3e-9;%[m]

%calc
%according to Thomas Ostler -> Zhifeng email (April.20th,2017),
%{
This seems like the usual overestimation of the Curie temperature using mean field.
I am not sure how clear it is in the paper but we corrected for the fact that the MF overestimate.
So, using your value of the Fe-Fe exchange as 4.5e-21 J, we use the expression for the Tc of an fcc lattice in the Heisenberg model:
3.18*Jat=kB*Tc, (Jat = atomistic exchange per interaction) which gives a Tc of 1037K. We then use the mean field expression:
JMF/3=kB*Tc, which gives a JMF (mean field exchange) of 4.29e-20 J.
%}
Tc_TM=3.18*Jat_TMTM/kb;
JMF_TMTM=kb*Tc_TM*3;
%similarly for RE
Tc_RE=3.18*Jat_RERE/kb;
JMF_RERE=kb*Tc_RE*3;

colT=size(T_,2);
totstep=round(runt/tstep);

szHext=size(Hextz_,2);
szstt=size(Jc_,2);
if calcFM
    mTMFM_=zeros(1,colT);mREFM_=zeros(1,colT);
elseif addSTT
    mREstt_=zeros(szHext,3,totstep+1);mTMstt_=zeros(szHext,3,totstep+1);
else
    mREHext_=zeros(szHext,3,totstep+1);mTMHext_=zeros(szHext,3,totstep+1);
end

for ctx=1:szx
    x=x_(ctx);
    q=1-x;%FM concentration
    Ms0=abs(q*Ms0_Fe-x*Ms0_Gd);
    Msperatom=abs(x*7.63-q*2.217);
    
    if x==0
        J0RERE=0;
        J0TMTM=JMF_TMTM;
        J0TMRE=0;
        J0RETM=0;
    elseif x==1
        J0RERE=JMF_RERE;%[J]
        J0TMTM=0;
        J0TMRE=0;
        J0RETM=0;
    else
        J0RERE=x*JMF_RERE;%[J]
        J0TMTM=q*JMF_TMTM;
        J0TMRE=x*z*JatTMRE;
        J0RETM=q*z*JatTMRE;
    end
    
    for cthext=1:szHext
        Hext=[0,0,Hextz_(cthext)];
        for ctstt=1:szstt
            Jc=Jc_(ctstt);
            for ctT=1:colT
                T=T_(ctT);
                bbeta=1/(kb*T);
                [mTTM,mTRE]=cweqn_wSTT(Hext,D,muRE,muTM,J0RERE,J0TMTM,J0TMRE,J0RETM,...
                    kb,T,x,q,mub,Msperatom,Ms0,ita,PFL,Jc,hbar,elev,tFL,alp,ip,lang_or_bri,JFe,JGd,addSTT,addSOT,thetaSHE);
                %mv:mTM; mk:mRE
                if calcFM
                    mTMFM_(ctT)=mTTM;mREFM_(ctT)=mTRE;
                else
                    if mREinitt
                        thet_init_TM=170*pi/180;
                    else
                        thet_init_TM=10*pi/180;
                    end
                    thet_init_RE=pi+thet_init_TM;
                    phi_init_TM=0;phi_init_RE=0;
                    mTM_init=[sin(thet_init_TM)*cos(phi_init_TM),sin(thet_init_TM)*sin(phi_init_TM),cos(thet_init_TM)];
                    mRE_init=[sin(thet_init_RE)*cos(phi_init_RE),sin(thet_init_RE)*sin(phi_init_RE),cos(thet_init_RE)];
                    
                    mTM_=zeros(totstep+1,3);mRE_=zeros(totstep+1,3);
                    mangle_=zeros(1,totstep+1);%angle between mTM and mRE
                    for ctrun=1:totstep
                        if ctrun==1
                            mTM=mTM_init;mRE=mRE_init;
                            mTM_(ctrun,:)=mTM_init;mRE_(ctrun,:)=mRE_init;
                        end
[Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA,...
    HRE_MFA,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,x,q,Msperatom,Ms0,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta);
                        ddt1_TM=LLB_integral(mTM,HTM_MFA,m0_TM,gam,Gam_parall_TM,Gam_perp_TM);
                        ddt1_RE=LLB_integral(mRE,HRE_MFA,m0_RE,gam,Gam_parall_RE,Gam_perp_RE);
                        
                        mTM=mTM_(ctrun,:)+ddt1_TM*tstep/2;
                        mRE=mRE_(ctrun,:)+ddt1_RE*tstep/2;
                        
[Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA,...
    HRE_MFA,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,x,q,Msperatom,Ms0,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta);
                        ddt2_TM=LLB_integral(mTM,HTM_MFA,m0_TM,gam,Gam_parall_TM,Gam_perp_TM);
                        ddt2_RE=LLB_integral(mRE,HRE_MFA,m0_RE,gam,Gam_parall_RE,Gam_perp_RE);
                        
                        mTM=mTM_(ctrun,:)+ddt2_TM*tstep/2;
                        mRE=mRE_(ctrun,:)+ddt2_RE*tstep/2;
                        
 [Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA,...
    HRE_MFA,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,x,q,Msperatom,Ms0,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta);
                        ddt3_TM=LLB_integral(mTM,HTM_MFA,m0_TM,gam,Gam_parall_TM,Gam_perp_TM);
                        ddt3_RE=LLB_integral(mRE,HRE_MFA,m0_RE,gam,Gam_parall_RE,Gam_perp_RE);
                        
                        mTM=mTM_(ctrun,:)+ddt3_TM*tstep;
                        mRE=mRE_(ctrun,:)+ddt3_RE*tstep;
                        
[Gam_parall_TM,Gam_parall_RE,Gam_perp_TM,Gam_perp_RE,HTM_MFA,...
    HRE_MFA,m0_TM,m0_RE]=heffcalc_LLB(D,muRE,muTM,mRE,mTM,mTTM,mTRE,...
    Hext,x,q,Msperatom,Ms0,ita,PFL,Jc,elev,tFL,J0RERE,J0TMTM,J0RETM,J0TMRE,alp,ip,bbeta);
                        ddt4_TM=LLB_integral(mTM,HTM_MFA,m0_TM,gam,Gam_parall_TM,Gam_perp_TM);
                        ddt4_RE=LLB_integral(mRE,HRE_MFA,m0_RE,gam,Gam_parall_RE,Gam_perp_RE);
                        
                        mTM=mTM_(ctrun,:)+tstep/6*(ddt1_TM+2*ddt2_TM+2*ddt3_TM+ddt4_TM);
                        mRE=mRE_(ctrun,:)+tstep/6*(ddt1_RE+2*ddt2_RE+2*ddt3_RE+ddt4_RE);
                        
                        mangle_(ctrun)=acos(dot(mTM,mRE)/(norm(mTM)*norm(mRE)));
                        mTM_(ctrun+1,:)=mTM;%no normalization because m isn't fixed length
                        mRE_(ctrun+1,:)=mRE;
                        
                        if (0)
                            mRE_(ctrun+1,:)=mRE_tmp/norm(mRE_tmp);
                        end
                    end
                    tt_=linspace(0,runt,totstep+1);
                end
            end
            if calcFM
                pplot2();
            elseif addSTT
                mREstt_(ctstt,:,:)=mRE_';mTMstt_(ctstt,:,:)=mTM_';
            else
                mREHext_(cthext,:,:)=mRE_';mTMHext_(cthext,:,:)=mTM_';
            end
        end
    end
end
if debbbug
    pplot();
elseif ppc==3
    save('final.mat')
end
toc