clear all
close all
global mod rhe
global R rho right result stressdrop pointlist %sliplength
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
generateMesh(model,'Hmax',500,'GeometricOrder','quadratic','Jiggle','on');
%% create model
c=1;
a=0;
f=0;
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);
% model.SolverOptions.ReportStatistics = 'on';
result = solvepde(model);
v = result.NodalSolution; %[m s^-1]
% elastic shear strain rate @ t=0
% esyx=result.XGradients;% [s^-1]
% esyz=result.YGradients;% [s^-1]
dt=3600*24*365*50;%[s]
% stressyx=result.XGradients*rhe.G*dt/2; %[Pa]
% stressyz=result.YGradients*rhe.G*dt/2; %[Pa]
% stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]

%% Pressure and Temperature
Tfield=result.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
Pfield=result.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
vsave=cell(2,1);
vsave{1,1}=v;
%% preset parameters
stressyx0=zeros(size(v));
stressyz0=zeros(size(v));
stress0=zeros(size(v));
esyxv=zeros(size(v));
esyzv=zeros(size(v));
stressdropyx=zeros(size(v));
stressdropyz=zeros(size(v));
faultstress0=zeros(mod.ny,1);
xq = linspace(0,0,mod.ny);
yq = linspace(0,mod.Y,mod.ny);
maxoverlimit=0;
stresslimit=byerlee(yq,'m')'*1e6;
loops=1890;
mov(loops) = struct('cdata',[],'colormap',[]);
savedata=cell(loops,3);
%%

for i=1:loops
    % total shear strain rate
    esyx=result.XGradients;% [s^-1]
    esyz=result.YGradients;% [s^-1]
    %% elastic stress accumulation(equal to the viscous stress)
    stressyx=stressyx0+(esyx-esyxv)*rhe.G*dt-stressdropyx; %[Pa] stress of previous time step; elastic stress accumulation; stress drop due to earthquake; -stressdropyx
    stressyz=stressyz0+(esyz-esyzv)*rhe.G*dt-stressdropyz; %[Pa]
    stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
    stressyx0=stressyx;
    stressyz0=stressyz;
    dstress=stress-stress0;
    stress0=stress;
    % k2: second invariant of strain rate tensor
    k2=((stress/1e6).^rhe.n./Fdisl).^2; %[s^-2]%
    etaeffdisl=(Fdisl).^(1/rhe.n).*1e6.*(k2).^((1-rhe.n)/(2*rhe.n)); %[Pa*s] effective viscosity for dislocation creep
    etaeff=etaeffdisl./2; % Assuming dislocation creep and diffusion creep has same effective viscosity;
    esyxv=stressyx./etaeff; % viscous shear strain rate
    %% viscous shear strain rate
    resultyx=createPDEResults(model,esyxv); % second order gradient of v_{v}
    resultyz=createPDEResults(model,esyzv); % second order gradient of v_{v}
    right=createPDEResults(model,resultyx.XGradients+resultyz.YGradients); % d2v_{v}/dx2+d2v_{v}/dz2
    specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
    %     model.SolverOptions.ReportStatistics = 'on';
    result = solvepde(model);
    %% total velocity
    v = result.NodalSolution; %[m s^-1]
    vsave{i+1,1}=v;
    %% Earthquake
    %  vdiff=vsave{i+1,1}-vsave{i,1};
    %  pdeplot(model,'xydata',dstress,'contour','on','levels',10)
    faultstress = interpolateSolution(createPDEResults(model,stressyx),xq,yq); % shear stress tau_{yx}
    %     dfaultstress = faultstress-faultstress0; % change of shear stress tau_{yx}
    %     faultstress0=faultstress;
    %     plot(yq,dfaultstress)
    %     title(['faultstress' num2str(i)])
    %     drawnow
    
    overlimit=faultstress-stresslimit;
    [overlimit,sliplength]=cutNegativeToZero(overlimit);
    plot(yq,overlimit,'g',yq,faultstress,'b',yq,stresslimit,'c');
    savedata{i,1}=overlimit;
    savedata{i,2}=faultstress;
    savedata{i,3}=stresslimit;
    maxoverlimit=max(overlimit);
    maxfaultstress=max(faultstress);
    minfaultstress=min(faultstress);
    %     ylim([0 maxfaultstress])
    %     xlim([0 maxfaultstress*0.85/1e4])
    title(['overlimit' num2str(i)])
    drawnow
    mov(i) = getframe(gcf);
    stressdrop=griddedInterpolant(yq,overlimit);
    %     sliplength=sliplength(1)*mod.Y/(mod.ny-1);
    slipregion0=(find(overlimit>0))*mod.Y/(mod.ny-1);
    slipregion=[slipregion0(1)-mod.Y/(mod.ny-1);slipregion0];
    disp([maxoverlimit maxfaultstress minfaultstress])
    %% solve for ƒ¢stress due to earthquake
    % find fault location
    cuts=diff(slipregion);
    cutsidx=find(cuts~=mod.Y/(mod.ny-1));
    regionidx1=cutsidx;
    regionidx2=cutsidx+1;
    if ~isempty(cutsidx)
        pointlist=[0;sort([slipregion(regionidx1);slipregion(regionidx2)-mod.Y/(mod.ny-1)]);slipregion(end);mod.Y];
    else
        pointlist=[0;slipregion(end);mod.Y];
    end
    nbs=length(pointlist)+2; % number of boundary segments
    %% create modeln
    modeln = createpde;
    gn=@npointrect;
    geometryFromEdges(modeln,gn);
    if length(pointlist)>3
        applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'u',0); %nonslip region at fault
        applyBoundaryCondition(modeln,'edge',nbs-3:nbs-1,'u',0);
        applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'g',@bcneumann);
    else
        applyBoundaryCondition(modeln,'edge',2:4,'u',0);
        applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann);
    end
    generateMesh(modeln,'Hmax',500,'GeometricOrder','quadratic','Jiggle','on');
    %     c=1;
    %     a=0;
    %     f=0;
    specifyCoefficients(modeln,'m',0,'d',0,'c',c,'a',a,'f',f);
    %     model6.SolverOptions.ReportStatistics = 'on';
    resultu = solvepde(modeln);
    stressdropyx=rhe.G*resultu.XGradients;
    stressdropyz=rhe.G*resultu.YGradients;
    stressdropyx=meshchange(model,modeln,stressdropyx);
    stressdropyz=meshchange(model,modeln,stressdropyz);
    %     if max(overlimit)
    %     end
    %     maxeta=max(etaeff(etaeff<inf));
    
    %   display([overlimit(100,1),maxeta,(i-1)*dt/(365*24*3600)]);
    
    %     if overlimit(100,1)>=0
    %       break
    %     end
    %     stressinterp=interpolateSolution(createPDEResults(model,stress),0,15000)
    %         maxstress=max(stress);
    %     if maxstress>3e8
    %         break
    %     end
    %     if stressinterp>byerlee(15000,'m')*1e6
    %         break
    %     end
end

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