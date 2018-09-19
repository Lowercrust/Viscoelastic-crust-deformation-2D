function tau=stresslowerlimit(z)
% input depth in meter
% output stress in Pa
global rho
tau=(rho*9.8*z*0.6);
end