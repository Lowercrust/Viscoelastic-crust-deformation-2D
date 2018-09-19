clear all
close all
global mod rhe
global R rho right result stressdroponfault pointlist overlimit region1 dt f1%sliplength
%% read data from file
filename='anorthite(wet).xlsx';
[~,~,model]=xlsread('model.xlsx');
[~,~,rheology]=xlsread(filename);
%% model parameter setting
mod=cell2table(model(:,2)');
mod.Properties.VariableNames=model(:,1)';
mod=table2struct(mod);
rhe=cell2table(rheology(:,2)');
rhe.Properties.VariableNames=rheology(:,1)';
rhe=table2struct(rhe);
R=8.3144598;%gas constant [J K^-1 mol^-1]
rho=2800;%density[kg m^-3]\
xq = linspace(0,0,mod.ny);
yq = linspace(0,mod.Y,mod.ny);
stresslimit=byerlee(yq,'m')'*1e6;
stresslowerlimit=(rho*9.8*yq*0.6)';
%% boundary and mesh
v0=mod.v0/(365*3600*24*1000);%[m/s]
% len=size(x1,2);
% if isempty(varargin)
%% rectangular model geometry
creepzone=50;
pointlist=[0;creepzone;mod.Y];
model = createpde;
gn=@npointrect;
geometryFromEdges(model,gn);
%% Boundary Conditions
applyBoundaryCondition(model,'edge',1,'g',@bcneumann1);
applyBoundaryCondition(model,'edge',4,'u',v0);
applyBoundaryCondition(model,'edge',2,'u',0);
generateMesh(model,'Hmax',500,'GeometricOrder','quadra','Jiggle','on','MesherVersion','R2013a');
%% create model
c=1;
a=0;
right=createPDEResults(model,zeros(size(model.Mesh.Nodes,2),1));
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
% model.SolverOptions.ReportStatiss = 'on';
% model.SolverOptions.AbsoluteTolerance=1e-60;
% model.SolverOptions.RelativeTolerance=1e-60;
% model.SolverOptions.ResidualTolerance=1e-60;
% model.SolverOptions.MinStep=1e-60;
% model.SolverOptions.MaxIterations=100;
% model.SolverOptions.ReportStatistics='on';
% elas shear strain rate @ t=0
% esyx=result.XGradients;% [s^-1]
% esyz=result.YGradients;% [s^-1]
dt=3600*24*365*1000;%[s]
dt0=dt;
elapsedtime=0;
% stressyx=result.XGradients*rhe.G*dt/2; %[Pa]
% stressyz=result.YGradients*rhe.G*dt/2; %[Pa]
% stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
%% Pressure and Temperature
Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
Pfield=model.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
%% preset parameters
stressyx0=zeros(size(model.Mesh.Nodes,2),1);
stressyz0=zeros(size(model.Mesh.Nodes,2),1);
stress0=zeros(size(model.Mesh.Nodes,2),1);
stressyx=zeros(size(model.Mesh.Nodes,2),1);
stressyz=zeros(size(model.Mesh.Nodes,2),1);
stress=zeros(size(model.Mesh.Nodes,2),1);
esyxv0=zeros(size(model.Mesh.Nodes,2),1);
esyzv0=zeros(size(model.Mesh.Nodes,2),1);
esyxv=zeros(size(model.Mesh.Nodes,2),1);
esyzv=zeros(size(model.Mesh.Nodes,2),1);
stressdropyx=zeros(size(model.Mesh.Nodes,2),1);
stressdropyz=zeros(size(model.Mesh.Nodes,2),1);
faultstress0=zeros(mod.ny,1);

maxoverlimit=0;

loops=10000;
earthquakes=100;
% mov(loops) = struct('cdata',[],'colormap',[]);
savedata=zeros(loops,7);
vsave=cell(loops,1);

% vsave{1,1}=v;
for j=1:earthquakes
    for i=1:loops
        %% stress Å®Å@effective viscosity and viscous shear strain rate
        [etaeff,esyxv,esyzv]=stress2visco(stress,Fdisl,stressyx,stressyz);
        %% calculate right hand side
        resultstressyx=createPDEResults(model,stressyx);
        resultstressyz=createPDEResults(model,stressyz);
        resultyx=createPDEResults(model,esyxv); % 2nd order gradient of v_{v}
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
        right=createPDEResults(model,right3); % d2v_{v}/dx2+d2v_{v}/dz2 -(resultstressyx.XGradients+resultstressyz.YGradients)./(rhe.G*dt)
        righttest=-interpolateSolution(right,region1.x,region1.y);
        %     disp(mean(righttest)/mean(f1))
        %% solve the equation for total shear strain rate
        result = solvepde(model);
        v = result.NodalSolution; %[m s^-1]
        % Total shear strain rate
        esyx=result.XGradients;% [s^-1]
        esyz=result.YGradients;% [s^-1]
        % verifiy
        a1=createPDEResults(model,esyx); %obtain 2nd order diff
        a2=createPDEResults(model,esyz); %obtain 2nd order diff
        left=a1.XGradients+a2.YGradients;
        remain=left-right.NodalSolution;
        %% stress accumulation from t_{0} to t_{0}+É¢t
        stressyx=stressyx+(esyx-esyxv)*rhe.G*dt; %[Pa] stress of previous time step; elas stress accumulation; stress drop due to earthquake; -stressdropyx
        stressyz=stressyz+(esyz-esyzv)*rhe.G*dt; %[Pa]
        stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
        %%
        %     specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
        %     model.SolverOptions.ReportStatiss = 'on';
        %     result = solvepde(model);
        %% Viscous shear strain rate
        % k2: second invariant of strain rate tensor
        
        %     pdeplot(model,'xydata',stress,'contour','on','levels',10)
        
        % plot(yq,faultstress,'b',yq,stresslimit,'c');
        % xlim([0 mod.X]);
        % ylim([0 max(stress)]);
        % %     %     maxoverlimit=max(overlimit);
        % title(['v' num2str(i)])
        % drawnow
        %     if max(overlimit)
        %     end
        %     maxeta=max(etaeff(etaeff<inf));
        
        %   display([overlimit(100,1),maxeta,(i-1)*dt/(365*24*3600)]);
        %     dt0=dt;
        elapsedtime=elapsedtime+dt;
        disp(elapsedtime/(365*24*3600))
        faultstress = interpolateSolution(createPDEResults(model,stressyx),xq,yq);
        plot(yq,faultstress,'r',yq,stresslimit,'c',yq,stresslowerlimit,'k')
        if max(faultstress)>0
            ylim([0,max(faultstress)]);
            xlim([0,max(faultstress)/(rho*9.8)]);
        end
        title(['interseismic' num2str(j) ',' num2str(i)])
        drawnow
%         stressyx0=stressyx;
%         stressyz0=stressyz;
%         stress0=stress;
        %         if max(stress>byerlee(15000,'m')*1e6)
        overlimit=faultstress-stresslimit;
        if ~isempty(overlimit(overlimit>0))
            break
        end
    end
    %% Eartbquakes
    acc=1;
    gn=@npointrect;
    xq1 = linspace(0,0,(mod.ny-1)*acc+1);
    yq1 = linspace(0,mod.Y,(mod.ny-1)*acc+1);
    faultstress = interpolateSolution(createPDEResults(model,stressyx),xq1,yq1);
    overlimit=faultstress-stresslowerlimit;
    % [overlimit,sliplength]=cutNegativeToZero(overlimit);
    stressdroponfault=griddedInterpolant(yq1,overlimit);
    depth0=find(diff(overlimit)<0);
    depth0=depth0(1)*50;
%     options = optimset('Display','iter');
    turn = fzero(@turningpoint,depth0);
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
    %     pdegplot(modeln, 'edgeLabels', 'on');
    applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
    applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann);
    % applyBoundaryCondition(modeln,'edge',2,'u',0);
    % applyBoundaryCondition(modeln,'edge',4,'g',@bcneumann);
    specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
    [u,mesh]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',10000,'ngen',inf); % coseismic slip
    resultn=createPDEResults(modeln,u);
    stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake
    stressdropyz=rhe.G*resultn.YGradients;
    stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2);
    stressdropyxc=meshchange(model,modeln,stressdropyx); % change modeln's mesh to model's mesh
    stressdropyzc=meshchange(model,modeln,stressdropyz);
    stressdropc=meshchange(model,modeln,stressdrop);
    stressyx=stressyx-stressdropyxc;
    stressyz=stressyz-stressdropyzc;
    faultstress = interpolateSolution(createPDEResults(model,stressyx),xq,yq); % shear stress tau_{yx}
    resultstressdroponfault = interpolateSolution(createPDEResults(model,stressdropc),xq,yq);
    resultstressdroponfaultyx = interpolateSolution(createPDEResults(model,stressdropyxc),xq,yq);
    resultstressdroponfaultyz = interpolateSolution(createPDEResults(model,stressdropyzc),xq,yq);
    plot(yq,faultstress,'r',yq,stresslimit,'c',yq,stresslowerlimit,'k')
    if max(faultstress)>0
        ylim([0,max(faultstress)]);
        xlim([0,max(faultstress)/(rho*9.8)]);
    end
    title(['After earthquake' num2str(j)])
    drawnow
    maxfaultstress=max(faultstress);
    minfaultstress=min(faultstress);
    dstress=stress-stress0;
    disp([maxfaultstress minfaultstress max(dstress) minetaeff turn dt/(365*24*3600) elapsedtime/(365*24*3600)])
end
%%
% stressdrop=(faultstress-cutstress);
% stressdrop=cutNegativeToZero(stressdrop);
% [vx2,~]=pdegrad(result.Mesh.Nodes,result.Mesh.Elements,esyxv);
% vx2=pdeprtni(result.Mesh.Nodes,result.Mesh.Elements,vx2);
% [~,vz2]=pdegrad(result.Mesh.Nodes,result.Mesh.Elements,esyzv);
% vz2=pdeprtni(result.Mesh.Nodes,result.Mesh.Elements,vz2);
% right=vx2+vz2;
%right = pdeInterpolant(result.Mesh.Nodes,result.Mesh.Elements,right); % create the interpolant
%specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
% x=0:dx:mod.X;
% y=0:dy:mod.Y;
% uxy = tri2grid(result.Mesh.Nodes,result.Mesh.Elements,right,x,y);