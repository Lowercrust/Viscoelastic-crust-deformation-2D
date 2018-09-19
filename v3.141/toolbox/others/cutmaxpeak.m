function [stress]=cutmaxpeak(model,stress)
global rho
% [~,I]=max(stress);
% [~,~,t]=meshToPet(model.Mesh);
% [~,col]=find(t(1:3,:)==I);
% ntl = pdeent(t,col);
% ntl = reshape(ntl,[length(col)*3,1]);
% meanstress=sum(stress(ntl))/length(ntl);
% meanstress1=sum(stress1(ntl))/length(ntl);
% stress(ntl)=meanstress;
% stress1(ntl)=meanstress1;
for i=1:length(model.Mesh.Nodes)
   if stress(i)>(model.Mesh.Nodes(2,i)*9.8*rho*0.6)
       stress(i)=model.Mesh.Nodes(2,i)*9.8*rho*0.6;
   end
end

% stress(I)=0;
% stress(I)=max(stress);
end