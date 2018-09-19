clear
close all
delete(gcp('nocreate'))
global mod rhe yieldstress ther
global R rho pointlist overlimit resultstressyx rockstrength faultstrength faultdepth %sliplength
formatOut = 'yymmddHHMMSS';
Simulationstarttime=datestr(now,formatOut);
savefolder=['E:\datasave\' Simulationstarttime ];
savealldatafolder=['E:\datasave\' Simulationstarttime '\alldata'];
% mkdir(savealldatafolder)
mkdir(savealldatafolder)
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
% applyBoundaryCondition(model,'edge',3,'u',v0);
% applyBoundaryCondition(model,'edge',[7,8],'u',0);
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
poolobj = parpool(6);
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
interseismiccount=0;
coseismiccount=0;
earthquakecount=1;
savefield=struct('tyx',{},'tyz',{},'result',{},'model',{},'Elapsedtime',{},'period',{});
for j=1:loops
    %% save data
%     Earthquake=struct('InterseismicStress',{},'Elapsedtime',{},'CoseismicStress',{},'CoseismicSlip',{});
    %(1) stressyx (2)stressyz (3)total shear strain rate yx (4) total shear
    %strain rate yz (5) interseismic or coseismic

 

    ii=ii+1;
    if interseismiccount==0 && coseismiccount==1
        earthquakecount=earthquakecount+1;
    end
    if isinterseismic(faultstress,faultnodesdepth)
        cut=0;
        coseismiccount=0;
        interseismiccount=interseismiccount+1;
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
        end
        tic
        [stress,stressyx,stressyz,elapsedtime,v,dt]=stressaccumulationplastic(stress,Fdisl,stressyx,stressyz,elapsedtime,model);
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
        title(['EarthquakeNo.' num2str(earthquakecount) 'Time' num2str(elapsedtime/(365*3600*24))])
        hold on
        scatter(faultnodesdepth,faultstress)
        hold off
        ylim([0,rho*9.8*faultdepth+100]);
        xlim([0,faultdepth*2+100]);
        drawnow
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=v;
        savefield(ii).model=model;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period='Interseismic';
    else
        interseismiccount=0;
        coseismiccount=coseismiccount+1;
        % update resultstressyx for specifiy boundary condition
        if exist('modeln','var') % if exist new model
            resultstressyx=createPDEResults(modeln,stressyx);
        else
            resultstressyx=createPDEResults(model,stressyx);
        end
        stresslowerlimit=(rho*9.8*faultnodesdepth*0.6)';
        overlimit=faultstress-stresslowerlimit;
        [pointlist,nbs]=findsection(overlimit,faultnodesdepth);
        faultdepth=pointlist(end-1);
        % interpolate old stress field to new meshgrid
        %% solve equation and update stress
        % input old stress field, out put new stress field and new model
        % and solution
        [stress,stressyx,stressyz,u,faultnodesdepth,faultstress,faultslip,modeln]=stressdrop(stressyx,stressyz,model);
        % new model become old model.
        model=modeln;
        % update resultstressyx
        %% plot fault stress
        [faultnodesdepth,faultstress,faultslip]=sortdata(faultnodesdepth,faultstress,faultslip);
        plot(faultnodesdepth,faultstress,'r',yq1,stresslimitplot,'c',yq1,stresslowerlimitplot,'k',faultnodesdepth,-faultslip*1e9,'g')
        hold on
        plot([faultdepth,faultdepth],[0,1e10],':');
        scatter(faultnodesdepth,faultstress)
        hold off
        ylim([0,rho*9.8*faultdepth]);
        xlim([0,faultdepth*2]);
        title(['EarthquakeNo.' num2str(earthquakecount) '@Time' num2str(elapsedtime/(365*3600*24)) '/EventNo.' num2str(coseismiccount)])
        drawnow
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=u;
        savefield(ii).model=modeln;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period=['coseismic No.' num2str(coseismiccount)];
    end
    if ii==100
        jj=jj+1;
        save([savefolder '\' num2str(jj*ii) '.mat'],'savefield');
        save([savealldatafolder  num2str(jj*ii) '.mat']);
        % clear the savefield and saver count
        savefield=struct('tyx',{},'tyz',{},'result',{},'model',{},'Elapsedtime',{},'period',{});
        ii=0;
    end
end
delete(poolobj)