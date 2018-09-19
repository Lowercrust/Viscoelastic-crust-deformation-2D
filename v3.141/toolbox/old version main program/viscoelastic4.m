clear
close all
delete(gcp('nocreate'))
global mod rhe yieldstress ther
global R rho pointlist overlimit resultstressyx rockstrength faultstrength faultdepth %sliplength
formatOut = 'yymmddHHMMSS';
Simulationstarttime=datestr(now,formatOut);
savefolder=['G:\datasave\' Simulationstarttime];
mkdir(savefolder)
%% read data from file
filename='anorthite(wet).xlsx';
[~,~,model]=xlsread('model.xlsx');
[~,~,thermal]=xlsread('thermal.xlsx');
[~,~,rheology]=xlsread(filename);
%% model parameter setting
mod=cell2table(model(:,2)');
mod.Properties.VariableNames=model(:,1)';
mod=table2struct(mod);
ther=cell2table(thermal(:,2)');
ther.Properties.VariableNames=thermal(:,1)';
ther=table2struct(ther);
rhe=cell2table(rheology(:,2)');
rhe.Properties.VariableNames=rheology(:,1)';
rhe=table2struct(rhe);
R=8.3144598;%gas constant [J K^-1 mol^-1]
rho=2800;%density[kg m^-3]\
yieldstress=1e8;
faultdepth=0; % initial faultdepth
faultstrength=1e6;
%% boundary and mesh
v0=mod.v0/(365*3600*24*1000);%[m/s]
% len=size(x1,2);
% if isempty(varargin)
%% rectangular model geometry (initial model without fault)
% model = createpde;
% R1 = [3,4,0,mod.X,mod.X,0,0,0,mod.Y,mod.Y]';
% geom = R1;
% % Names for the geometric object
% ns = (char('R1'))';
% % Set formula
% sf = 'R1';
% % Create geometry
% g = decsg(geom,sf,ns);
g=@rectwithsubdomain;
model = createpde;
geometryFromEdges(model,g);
%% Boundary Conditions
applyBoundaryCondition(model,'edge',1,'u',v0);
applyBoundaryCondition(model,'edge',[7,8],'u',0);
%% Coeffecient
% faultstress0=zeros(mod.ny,1);
xq1 = linspace(0,0,mod.Y-1+1); % 1m grid
yq1 = linspace(0,mod.Y,mod.Y+1);
maxoverlimit=0;
rockstrength=1e6; % pa
stresslimitplot=(rho*9.8*yq1*0.6)'+rockstrength;
stresslowerlimitplot=(rho*9.8*yq1*0.6)';
loops=100000;
Events=struct;
poolobj = parpool(16);
overupperlimit=0;
dt=3600*24*365*1000;%[s]
dt0=dt;
elapsedtime=0;
% faultstress=zeros(size(yq1,2),1);
% faultslip=faultstress;
%% generateMesh for interseismic model (initial)
generateMesh(model,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
%% variables
nodesize=size(model.Mesh.Nodes,2);
% old stress field
stressyx0=zeros(nodesize,1);
stressyz0=zeros(nodesize,1);
stress0=zeros(nodesize,1);
% new stress field
stressyx=zeros(nodesize,1);
stressyz=zeros(nodesize,1);
stress=zeros(nodesize,1);
% old shear strain rate
esyxv0=zeros(nodesize,1);
esyzv0=zeros(nodesize,1);
% new shear strain rate
esyxv=zeros(nodesize,1);
esyzv=zeros(nodesize,1);
% stress drop
stressdropyx=zeros(nodesize,1);
stressdropyz=zeros(nodesize,1);
resultstressyx=createPDEResults(model,stressyx);
%% Pressure and Temperature (1-D linear field)
Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
Pfield=model.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
%% calculate faultstress and overupperlimit
faultnodes=find(model.Mesh.Nodes(1,:)==0); % find nodes on the fault
faultnodesdepth=model.Mesh.Nodes(2,faultnodes);
faultstress=stressyx(faultnodes);
% sort data
[faultnodesdepth,faultstress]=sortdata(faultnodesdepth,faultstress);
ii=0;
jj=0;
for j=1:loops
    %% save data
%     Earthquake=struct('InterseismicStress',{},'Elapsedtime',{},'CoseismicStress',{},'CoseismicSlip',{});
    %(1) stressyx (2)stressyz (3)total shear strain rate yx (4) total shear
    %strain rate yz (5) interseismic or coseismic
    savefield=struct('tyx',{},'tyz',{},'result',{},'model',{},'Elapsedtime',{},'period',{});
    ii=ii+1;
    if isinterseismic(faultstress,faultnodesdepth)
        
        % updatestress
        if exist('modeln','var')
            %reconstruct model
            model = createpde;
            model.Geometry=modeln.Geometry;
            applyBoundaryCondition(model,'edge',1,'u',v0);
            applyBoundaryCondition(model,'edge',[7,8],'u',0);
            model.Mesh=modeln.Mesh;
            Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
            Pfield=model.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
            Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
%             applyBoundaryCondition(model,'edge',2,'u',v0);
%             applyBoundaryCondition(model,'edge',4,'u',0);
        end
        tic
        [stress,stressyx,stressyz,elapsedtime,v]=stressaccumulation(stress,Fdisl,stressyx,stressyz,elapsedtime,model);
        timeforcal=toc;
        disp(['timeforcal.:' num2str(timeforcal) '/elapsedtime:' num2str(elapsedtime/(365*3600*24))])
        %% plot stress
        faultnodes=find(model.Mesh.Nodes(1,:)==0); % find nodes on the fault
        faultnodesdepth=model.Mesh.Nodes(2,faultnodes);
        faultstress=stressyx(faultnodes);
        % sort data
        [faultnodesdepth,faultstress]=sortdata(faultnodesdepth,faultstress);
        plot(faultnodesdepth,faultstress,'r',yq1,stresslimitplot,'c',yq1,stresslowerlimitplot,'k')
        ylim([0,max(faultstress)*2]);
        xlim([0,max(faultstress)*2/(rho*9.8)]);
        drawnow
%         Earthquake(ii).InterseismicStress=faultstress;
%         Earthquake(ii).Elapsedtime=elapsedtime/(365*24*3600);
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=v;
        savefield(ii).model=model;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period='Interseismic';
    else
        %% create new model with fault region
        gn=@rectwithsubdomain;
        if exist('modeln','var')
            resultstressyx=createPDEResults(modeln,stressyx);
        else
            resultstressyx=createPDEResults(model,stressyx);
        end
        stresslowerlimit=(rho*9.8*faultnodesdepth*0.6)';
        overlimit=faultstress-stresslowerlimit;
        [pointlist,nbs]=findsection(overlimit,faultnodesdepth);
        faultdepth=pointlist(end-1);
        %% create new model
        modeln = createpde;
        geometryFromEdges(modeln,gn);
        generateMesh(modeln,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
        % refine mesh for fault region
        modeln=refinemeshface(modeln,2,0.1);
        %% apply boundary condition
        %boundary condition on fault
%         applyBoundaryCondition(modeln,'edge',8,'g',@bcneumann,'Vectorized','off');
%         applyBoundaryCondition(modeln,'edge',[3,4,9,10],'u',0);
        applyBoundaryCondition(modeln,'edge',7,'g',@bcneumann,'Vectorized','off');
        applyBoundaryCondition(modeln,'edge',[1,2,8],'u',0);
        specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
        % interpolate old stress field to new meshgrid
        tic
        stressyx=meshchange(modeln,model,stressyx);
        toc
        stressyz=meshchange(modeln,model,stressyz);
        toc
        model=modeln;
        %% solve equation and update stress
        [stress,stressyx,stressyz,u]=stressdrop(stressyx,stressyz,modeln);
        % update resultstressyx
        
        %% plot fault stress
        faultnodes=find(modeln.Mesh.Nodes(1,:)==0); % find nodes on the fault
        faultnodesdepth=modeln.Mesh.Nodes(2,faultnodes);
        faultstress=stressyx(faultnodes);
        faultslip=u(faultnodes);
%         fault=sortrows([faultnodesdepth',faultstress,faultslip]);
        [faultnodesdepth,faultstress,faultslip]=sortdata(faultnodesdepth,faultstress,faultslip);
        plot(faultnodesdepth,faultstress,'r',yq1,stresslimitplot,'c',yq1,stresslowerlimitplot,'k',faultnodesdepth,-faultslip*1e9,'g')
        hold on
        plot([faultdepth,faultdepth],[0,1e10],':');
        hold off
        ylim([0,max(faultstress)*2]);
        xlim([0,max(faultstress)*2/(rho*9.8)]);
        drawnow
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=u;
        savefield(ii).model=modeln;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period='coseismic';
    end
    if ii==100
        jj=jj+1;
        save([savefolder '\' num2str(jj*ii) '.mat'],'savefield');
        ii=0;
    end
end
delete(poolobj)