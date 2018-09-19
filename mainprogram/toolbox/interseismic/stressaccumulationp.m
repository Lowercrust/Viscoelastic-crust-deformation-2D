function [stressyx,stressyz,etaeff,esyxv,esyzv,elapsedtime,v,dt,varargout]=stressaccumulationp(Fdisl,stressyx,stressyz,elapsedtime,model,femodel)
global modc fi stresslimitall
%% stress ??�??�??�@effective viscosity and viscous shear strain rate
[etaeff,esyxv,esyzv]=stress2visco(Fdisl,model.Mesh,stressyx,stressyz);
% resultstressyx=createPDEResults(model,stressyx);
% resultstressyz=createPDEResults(model,stressyz);
% resultyx=createPDEResults(model,esyxv); % 2nd order gradient of v_{v}
% resultyz=createPDEResults(model,esyzv); % 2nd order gradient of v_{v}

minetaeff=min(etaeff);
% disp(minetaeff)
if  isinf(minetaeff)
    dt=3600*24*365*50;%[s] first time step
else
    dt=minetaeff*0.1/modc.rhe{1}.G;
    if dt>3600*24*365*50
        dt=3600*24*365*50;
    end
end
% 
%% solve the equation for total shear strain rate
right=createright(model.Mesh,dt,'tyx',stressyx,'tyz',stressyz,'esyxv',esyxv,'esyzv',esyzv);
femodel=createfemodel(model,'update',right,femodel); % in each time step only right hand side changes
% K=distributed(femodel.K);
% F=distributed(femodel.F);
% tic;KF=K\F;toc
% KF=gather(KF);
v= femodel.B*(femodel.K\femodel.F) +femodel.ud;
% v=solvelinearPDE(model,femodel);
% [p,~,t] = meshToPet(model.Mesh);
p=model.Mesh.Nodes;
t=model.Mesh.Elements;
[esyx,esyz]=trigrad(p,t,v,'form','node');
%% stress accumulation from t_{0} to t_{0}+dt
if size(modc.rhe,1)==2
    for i=1:2
        stressyx(fi{i},1)=stressyx(fi{i},1)+(esyx(fi{i},1)-esyxv(fi{i},1))*modc.rhe{i}.G*dt;
        stressyz(fi{i},1)=stressyz(fi{i},1)+(esyz(fi{i},1)-esyzv(fi{i},1))*modc.rhe{i}.G*dt;
    end
else
    stressyx=stressyx+(esyx-esyxv)*rhe.G*dt; %[Pa] stress of previous time step; elas stress accumulation; stress drop due to earthquake; -stressdropyx
    stressyz=stressyz+(esyz-esyzv)*rhe.G*dt; %[Pa]
end
% if modc.considerplastisity
%     dstress=stressyx-stresslimitall;
%     dstress(dstress<0)=0;
%     stressyx=stressyx-dstress;
%     varargout{1}=dstress;
% end
elapsedtime=elapsedtime+dt;
end