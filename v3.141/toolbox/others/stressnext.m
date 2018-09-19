function [stressyx,stressyz]=stressnext(stressyx0,stressyz0,model)
global mod rhe yieldstress ther
global R rho right result resultn pointlist overlimit region1 resultstressyx  %sliplength
resultstressyx=createPDEResults(model,stressyx0);
resultstressyz=createPDEResults(model,stressyz0);
stress0=sqrt(stressyx0.^2+stressyz0.^2);
parfor ii=1:size(yq1,2)
    faultstress(ii) = interpolateSolution(resultstressyx1,xq1(ii),yq1(ii));
end
overupperlimit=faultstress-stresslimit;
if isempty(overupperlimit(overupperlimit>0))% interseismic
    %% stress Å®Å@effective viscosity and viscous shear strain rate
    [etaeff,esyxv,esyzv]=stress2visco(stress0,Fdisl,stressyx0,stressyz0);
    %% calculate right hand side
    resultstressyx=createPDEResults(model,stressyx);
    resultstressyz=createPDEResults(model,stressyz);
    resultyx=createPDEResults(model,esyxv); % 2nd order gradient of v_{v} viscosflow
    resultyz=createPDEResults(model,esyzv); % 2nd order gradient of v_{v}
    right1=(resultstressyx.XGradients+resultstressyz.YGradients)/(rhe.G*dt);
    minetaeff=min(etaeff);
    if  isinf(minetaeff)
        dt=3600*24*365*100;%[s]
    else
        dt=minetaeff*0.1/rhe.G;
    end
    right2=resultyx.XGradients+resultyz.YGradients;
    right3=right2-right1;
    right=createPDEResults(model,right3);
else % coseismic
    
end
end
