function [T,varargout]=geothermal(varargin)
global modc
% Qr=0.030;% [W/m^2]
Q0=0.065; % surface heat flow
if length(varargin)==1 % determine D 
    z=modc.mod.Y;
    D=varargin{1};
elseif length(varargin)==2 % calculate T
    z=varargin{1};
    D=varargin{2}; %[m]
end
h=modc.mod.CT;
% h=45000;
k1=modc.ther{1}.k; %W/(m*k)
k2=modc.ther{2}.k;
Qr=(Q0*exp(-h/D)/k1)/(1/k2-1/k1+exp(-h/D)/k1);
A=(Q0-Qr)/D;
varargout{1}=A;
% crustzi=find(z<modc.mod.CT);
% mantlezi=find(z>modc.mod.CT);
T=zeros(length(z),1);
T(z<=modc.mod.CT,1)=Qr.*z(z<=modc.mod.CT)./k1+D.^2.*A.*(1-exp(-z(z<=modc.mod.CT)./D))./k1;

Th=Qr.*h./k1+D.^2.*A.*(1-exp(-h./D))./k1;
% dThdz=Qr/k1+D*A*exp(-h/D)/k1;
% T(z>modc.mod.CT,1)=Th+z(z>modc.mod.CT)*
T(z>modc.mod.CT,1)=Th+(z(z>modc.mod.CT)-h).*Qr./k2;
T=T'+modc.mod.Ts;
end