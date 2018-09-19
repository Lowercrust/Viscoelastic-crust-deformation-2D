clear all
close all
global mod rhe
global R rho right result stressdroponfault pointlist overlimit region1 f1%sliplength
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
%% boundary and mesh
v0=mod.v0/(365*3600*24*1000);%[m/s]

% len=size(x1,2);
% if isempty(varargin)
%% rectangular model geometry
model = createpde;
R1 = [3,4,0,mod.X,mod.X,0,0,0,mod.Y,mod.Y]';
geom = R1;
% Names for the geometric object
ns = (char('R1'))';
% Set formula
sf = 'R1';
% Create geometry
g = decsg(geom,sf,ns);
geometryFromEdges(model,g);
%% Boundary Conditions
applyBoundaryCondition(model,'edge',2,'u',v0);
applyBoundaryCondition(model,'edge',4,'u',0);
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
dt=3600*24*365*100;%[s]
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
xq = linspace(0,0,mod.ny);
yq = linspace(0,mod.Y,mod.ny);
maxoverlimit=0;
stresslimit=byerlee(yq,'m')'*1e6;
loops=20000;
% mov(loops) = struct('cdata',[],'colormap',[]);
savedata=zeros(loops,7);
vsave=cell(loops,1);
% vsave{1,1}=v;
%%
result = solvepde(model);
v = result.NodalSolution; %[m s^-1]
%% Total shear strain rate
esyx=result.XGradients;% [s^-1]
esyz=result.YGradients;% [s^-1]
a1=createPDEResults(model,esyx);
a2=createPDEResults(model,esyz);
left=a1.XGradients+a2.YGradients;
remain=left-right.NodalSolution;
%% Elastic stress accumulation (equal to the viscous stress)
for j=1:1
    for i=1:loops
        stressyx=stressyx0+(esyx-esyxv)*rhe.G*dt; %[Pa] stress of previous time step; elas stress accumulation; stress drop due to earthquake; -stressdropyx
        stressyz=stressyz0+(esyz-esyzv)*rhe.G*dt; %[Pa]
        resultstressyx=createPDEResults(model,stressyx0);
        resultstres
        syz=createPDEResults(model,stressyz0);
        %%
        resultyx=createPDEResults(model,esyxv); % second order gradient of v_{v}
        resultyz=createPDEResults(model,esyzv); % second order gradient of v_{v}
        %     right1=(resultstressyx.XGradients+resultstressyz.YGradients)/(rhe.G*dt);
        %     right2=resultyx.XGradients+resultyz.YGradients;
        %     right3=right2-right1;
        %     right=createPDEResults(model,right2); % d2v_{v}/dx2+d2v_{v}/dz2 -(resultstressyx.XGradients+resultstressyz.YGradients)./(rhe.G*dt)
        %     specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
        %     model.SolverOptions.ReportStatiss = 'on';
        %     result = solvepde(model);
        %% Viscous shear strain rate
        % k2: second invariant of strain rate tensor
        stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
        [etaeff,esyxv,esyzv]=stress2visco(stress,Fdisl);
        minetaeff=min(etaeff);
        faultstress = interpolateSolution(createPDEResults(model,stressyx),xq,yq); % shear stress tau_{yx}
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
        if  isinf(minetaeff)
            dt=3600*24*365*100;%[s]
        else
            dt=minetaeff*0.1/rhe.G;
        end
        elapsedtime=elapsedtime+dt;
        maxfaultstress=max(faultstress);
        minfaultstress=min(faultstress);
        dstress=stress-stress0;
        disp([maxfaultstress minfaultstress max(dstress) minetaeff dt/(365*24*3600) elapsedtime/(365*24*3600)])
        stressyx0=stressyx;
        stressyz0=stressyz;
        stress0=stress;
        if max(stress>byerlee(15000,'m')*1e6)
            break
        end
    end
    %% Earthquake
    xq1 = linspace(0,0,mod.ny*100);
    yq1 = linspace(0,mod.Y,mod.ny*100);
    faultstress = interpolateSolution(createPDEResults(model,stressyx),xq1,yq1);
    stresslimit=byerlee(yq1,'m')'*1e6;
    overlimit=faultstress-stresslimit;
    [overlimit,sliplength]=cutNegativeToZero(overlimit);
    stressdroponfault=griddedInterpolant(yq1,overlimit);
    slipregion0=(find(overlimit>0))*mod.Y/(mod.ny-1);
    slipregion=[slipregion0(1)-mod.Y/(mod.ny-1);slipregion0];
    cuts=diff(slipregion);
    cutsidx=find(cuts~=mod.Y/(mod.ny-1));
    regionidx1=cutsidx;
    regionidx2=cutsidx+1;
    pointlist=[0;slipregion(end);mod.Y];
    nbs=length(pointlist)+2;
    modeln = createpde;
    gn=@npointrect;
    geometryFromEdges(modeln,gn);
    applyBoundaryCondition(modeln,'edge',2:4,'u',0);
    applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann);
    generateMesh(modeln,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
    specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
    resultu = solvepde(modeln);
    stressdropyx=rhe.G*resultu.XGradients;
    stressdropyz=rhe.G*resultu.YGradients;
    stressdropyx=meshchange(model,modeln,stressdropyx);
    stressdropyz=meshchange(model,modeln,stressdropyz);
    %% Post earthquake
    figure('Position',[0 100 800 600])
    for j=1:1
        stressyx=stressyx-stressdropyx;
        stressyz=stressyz-stressdropyz;
        if j==1
            stressyxi=stressyx-stressdropyx;
            stressyzi=stressyz-stressdropyz;
        end
        stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
        [etaeff,esyxv,esyzv]=stress2visco(stress,Fdisl,stressyx,stressyz);
        minetaeff=min(etaeff);
        dt=minetaeff*0.1/rhe.G;
        resultstressafteryx=createPDEResults(model,stressyx);
        resultstressafteryz=createPDEResults(model,stressyz);
        resultafteryx=createPDEResults(model,esyxv); % second order gradient of v_{v}
        resultafteryz=createPDEResults(model,esyzv); % second order gradient of v_{v}
        right=resultafteryx.XGradients+resultafteryz.YGradients-(resultstressafteryx.XGradients+resultstressafteryz.YGradients)/(rhe.G*dt);
        right=createPDEResults(model,right);
        result = solvepde(model);
        v = result.NodalSolution; %[m s^-1]
        esyx=result.XGradients;% [s^-1]
        esyz=result.YGradients;% [s^-1]
        stressyx=stressyx+(rhe.G*dt)*(esyx-esyxv);
        stressyz=stressyz+(rhe.G*dt)*(esyz-esyzv);
        subplot(2,2,1)       
        pdeplot(model,'xydata',stressyx,'contour','on','levels',10)
        title(['stressyx' num2str(j)])
        subplot(2,2,2)       
        pdeplot(model,'xydata',stressyxi,'contour','on','levels',10)
        title('stressyxbefore')
        subplot(2,2,3)       
        pdeplot(model,'xydata',stressyx-stressyxi,'contour','on','levels',10)
        title('diff')
        subplot(2,2,4)
        pdeplot(model,'xydata',v,'contour','on','levels',10)
        title('v')
        drawnow
    end
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