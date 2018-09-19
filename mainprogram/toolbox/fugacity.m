% caculation of water fugacity with van der Waals equation
function f=fugacity(P,T)
% caculate the fugacity of water using Van der Waals equation
% input unit [Pa] [K]
% output fugacity in the unit of [MPa]
global modc
a=5.536*0.1;
b=0.03049*0.001;
P0=0;
R=modc.R;
f=(exp(b.*(P-P0)./(R.*T))./((a.*P+R.^2*T.^2)./(a.*P0+R.^2*T.^2)).*P)./(1e6);
end