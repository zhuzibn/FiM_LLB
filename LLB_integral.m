
function ddt=LLB_integral(m,Heff,m0,gam,Gam_parall,Gam_perp)
ddt=gam*cross(m,Heff)-Gam_parall*(1-dot(m,m0)/(norm(m)^2))*m...
    -Gam_perp*cross(m,cross(m,m0))/(norm(m)^2);
% ddt=gam*cross(m,Heff)-Gam_parall*(1-dot(m,m0))*m...
%     -Gam_perp*cross(m,cross(m,m0));
end