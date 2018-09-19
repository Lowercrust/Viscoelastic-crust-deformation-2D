clear
close all
load('G:\result10241')
%% Earthquake
% stress drop on fault
acc=1;
gn=@npointrect;
xq1 = linspace(0,0,(mod.ny-1)*acc+1);
yq1 = linspace(0,mod.Y,(mod.ny-1)*acc+1);
faultstress = interpolateSolution(createPDEResults(model,stressyx),xq1,yq1);
stresslimit=byerlee(yq1,'m')'*1e6;
overlimit=faultstress-stresslimit;
% [overlimit,sliplength]=cutNegativeToZero(overlimit);
stressdroponfault=griddedInterpolant(yq1,overlimit);
depth0=0;
options = optimset('Display','iter');
turn = fzero(@turningpoint,depth0,options);
overlimit(overlimit<0)=0;
stressdroponfault=griddedInterpolant(yq1,overlimit);
% % slip region
% slipregion=(find(overlimit>0))*mod.Y/((mod.ny-1)*acc+1)-mod.Y/((mod.ny-1)*acc+1);
% % slipregion=[slipregion0(1)-mod.Y/((mod.ny-1)*acc+1);slipregion0];
% cuts=diff(slipregion);
% cutsidx=find(cuts~=mod.Y/((mod.ny-1)*acc+1));
% regionidx1=cutsidx;
% regionidx2=cutsidx+1;
pointlist=[0;turn;mod.Y];
% nbs=length(pointlist)+2;
% create model
modeln = createpde;
geometryFromEdges(modeln,gn);
% geometryFromEdges(modeln,g);
pdegplot(modeln, 'edgeLabels', 'on');
applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann);
% applyBoundaryCondition(modeln,'edge',2,'u',0);
% applyBoundaryCondition(modeln,'edge',4,'g',@bcneumann);
specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
[u,mesh]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',50000,'ngen',inf); % coseismic slip
resultn=createPDEResults(modeln,u);
stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake
stressdropyz=rhe.G*resultn.YGradients;
stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2);
stressdropyxc=meshchange(model,modeln,stressdropyx); % change modeln's mesh to model's mesh
stressdropyzc=meshchange(model,modeln,stressdropyz);
stressdropc=meshchange(model,modeln,stressdrop);

resultstressdroponfault = interpolateSolution(createPDEResults(model,stressdropc),xq,yq);
resultstressdroponfaultyx = interpolateSolution(createPDEResults(model,stressdropyxc),xq,yq);
resultstressdroponfaultyz = interpolateSolution(createPDEResults(model,stressdropyzc),xq,yq);
% stressdropyxc(stressdropyxc<0)=0;
% stressdropyzc(stressdropyzc<0)=0;
% stressdropc(stressdropc<0)=0;
% figure('Position',[0 100 800 600])
% for j=1:1
%     stressafteryx=stressyx-stressdropyxc;
%     stressafteryz=stressyz-stressdropyzc;
%     if j==1
%         stressafteryxi=stressyx-stressdropyxc;
%         stressafteryzi=stressyz-stressdropyzc;
%         subplot(2,2,2)
%         pdeplot(model,'xydata',stressafteryxi,'contour','on','levels',10)
%         title('stressyxbefore')
%     end
%     stressafter=sqrt(stressafteryx.^2+stressafteryz.^2); %[Pa]
%     [etaeff,esyxv,esyzv]=stress2visco(stressafter,Fdisl,stressafteryx,stressafteryz);
%     minetaeff=min(etaeff);
%     dt=minetaeff*0.1/rhe.G;
%     resultstressafteryx=createPDEResults(model,stressafteryx);
%     resultstressafteryz=createPDEResults(model,stressafteryz);
%     resultafteryx=createPDEResults(model,esyxv); % second order gradient of v_{v}
%     resultafteryz=createPDEResults(model,esyzv); % second order gradient of v_{v}
%     right=resultafteryx.XGradients+resultafteryz.YGradients; %-(resultstressafteryx.XGradients+resultstressafteryz.YGradients)/(rhe.G*dt)
%     right=createPDEResults(model,right);
%     result = solvepde(model);
%     v = result.NodalSolution; %[m s^-1]
%     esyx=result.XGradients;% [s^-1]
%     esyz=result.YGradients;% [s^-1]
%     resultesafteryx=createPDEResults(model,esyx);
%     resultesafteryz=createPDEResults(model,esyz);
%     left=resultesafteryx.XGradients+resultesafteryz.YGradients;
%     stressafteryx=stressafteryx+(rhe.G*dt)*(esyx-esyxv);
%     stressafteryz=stressafteryz+(rhe.G*dt)*(esyz-esyzv);
%     subplot(2,2,1)
%     pdeplot(model,'xydata',stressafteryx,'contour','on','levels',10)
%     title(['stressyx' num2str(j)])
%     
%     subplot(2,2,3)
%     pdeplot(model,'xydata',stressafteryx-stressafteryxi,'contour','on','levels',10)
%     title('diff')
%     subplot(2,2,4)
%     pdeplot(model,'xydata',v,'contour','on','levels',10)
%     title('v')
%     drawnow
% end