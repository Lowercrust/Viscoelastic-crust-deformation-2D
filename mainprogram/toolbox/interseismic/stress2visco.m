function [etaeff,esyxv,esyzv]=stress2visco(Fdisl,mesh,stressyx,stressyz)
global modc fi
stress=hypot(stressyx,stressyz);
if length(modc.rhe)==2
    crhe=modc.rhe{1};
    mrhe=modc.rhe{2};
    p=mesh.Nodes;
    p1i=fi{1};
    p2i=fi{2};
    k2=zeros(length(p),1);
    etaeffdisl=zeros(length(p),1);
    k2(p1i,1)=((stress(p1i,1)/1e6).^crhe.n./Fdisl(p1i,1)).^2; %[s^-2]%
    k2(p2i,1)=((stress(p2i,1)/1e6).^mrhe.n./Fdisl(p2i,1)).^2; %[s^-2]%
    etaeffdisl(p1i)=(Fdisl(p1i,1)).^(1/crhe.n).*1e6.*(k2(p1i,1)).^((1-crhe.n)/(2*crhe.n));
    etaeffdisl(p2i)=(Fdisl(p2i,1)).^(1/mrhe.n).*1e6.*(k2(p2i,1)).^((1-mrhe.n)/(2*mrhe.n));
elseif length(modc.rhe)==1
    crhe=modc.rhe{1};
    k2=((stress/1e6).^crhe.n./Fdisl).^2; %[s^-2]%
    etaeffdisl=(Fdisl).^(1/crhe.n).*1e6.*(k2).^((1-crhe.n)/(2*crhe.n)); %[Pa*s] effective viscosity for dislocation creep
end
etaeff=etaeffdisl./2; % Assuming dislocation creep and diffusion creep has same effective viscosity;
esyxv=2*stressyx./etaeff; % viscous shear strain rate yx component
esyzv=2*stressyz./etaeff; % viscous shear strain rate yz component
end