function f = shearheating(mesh,pres,Temp,esyxv,esyzv)
% calculate shear heating
% input arg.: mesh, temperature field, viscos shear strain rate
global modc fi
% f = zeros(N,nr); % Allocate f
p = mesh.Nodes;
t = mesh.Elements;
if length(fi)==2
    f=zeros(length(p),1);
    Fpt=FPT(pres,Temp,mesh,'disl')';
    for i=1:2
        n=modc.rhe{i}.n;
        f(fi{i},1)=0.25^(1/(2*n)).*(Fpt(fi{i},1)).^(1/n).*(esyxv(fi{i},1).^2+esyzv(fi{i},1).^2).^((1+n)/(2*n));
    end
else
    n = modc.rhe{1}.n;
    % p = modc.ther{1}.rho*9.8*mesh.Nodes(2,:)+modc.atm;
    
    f= 0.25^(1/(2*n)).*(FPT(pres,Temp,'disl')').^(1/n).*(esyxv.^2+esyzv.^2).^((1+n)/(2*n));
    % Interpolate from node data to triangle midpoint data
end
f= pdeintrpgpu(p,t,f);
end