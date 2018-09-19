function tau=stresslimit(y)
% input depth in meter
% output stress in Pa
global rho rockstrength resultstressyx
faultstrength=(rho*9.8*y*0.6)+rockstrength;
faultstress = interpolateSolution(resultstressyx,0,y);
tau=faultstress-faultstrength;
end