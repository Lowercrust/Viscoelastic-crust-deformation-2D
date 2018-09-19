clear
close all
addpath(genpath('toolbox'))
% delete(gcp('nocreate'))
% poolobj = parpool(4);
global mod rhe yieldstress ther atm
global R rho pointlist overlimit resultstressyx rockstrength faultstrength faultdepth faultwidth stressyxi%sliplength
formatOut = 'yymmddHHMMSS';
Simulationstarttime=datestr(now,formatOut);
ver=3.14;
addpath(genpath('toolbox'))
savefolder=['..' filesep 'datasave' filesep 'v' num2str(ver) Simulationstarttime];
mkdir(savefolder)
%% read data from file
[~,~,model]=xlsread(lsinputfile('Model Properties'));
[~,~,thermal]=xlsread(lsinputfile('Thermal Properties'));
[~,~,anorthite]=xlsread(lsinputfile('Rheology'));
%% model parameter setting
mod=cell2table(model(:,2)');
mod.Properties.VariableNames=model(:,1)';
mod=table2struct(mod);
% mod.v0=15;
ther=cell2table(thermal(:,2)');
ther.Properties.VariableNames=thermal(:,1)';
ther=table2struct(ther);
rhe=cell2table(anorthite(:,2)');
rhe.Properties.VariableNames=anorthite(:,1)';
rhe=table2struct(rhe);
R=8.3144598;%gas constant [J K^-1 mol^-1]
rho=2800;%density[kg m^-3]\
yieldstress=1e8;
faultdepth=25000; % initial faultdepth
faultstrength=5e6;
faultwidth=0.5;
atm=101325; % [pa] 
v0=mod.v0/(365*3600*24*1000);%[m/s] boundary velocity
% load mesh
load(lsmeshfile);
%% out put basic information for the data set
createmodelinfo(msh,mod,rhe,Simulationstarttime)
%% interseismic model
g=@rectwithsubdomain;
model = createpde;
geometryFromEdges(model,g);
% Boundary Conditions
applyBoundaryCondition(model,'edge',1,'u',v0);% constant velocity boundary
applyBoundaryCondition(model,'edge',[7,8],'u',0);
% adjust faultdepth
faultnodesindex=find(msh.Nodes(1,:)==0);
faultnodesy=msh.Nodes(2,faultnodesindex);
faultdepth=min(faultnodesy(faultnodesy>0));
% boundary condition matrix
%% coseismic model
gn=@fivepointrect;
modeln = createpde;
geometryFromEdges(modeln,gn);
% Boundary Conditions
applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','off');
applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
% refine mesh near surface 
%% temperature model
gt=@rectwithsubdomain;
modelt = createpde;
geometryFromEdges(modelt,gt);
% Boundary Conditions
thermalgradient = @(region,state)0.025*region.y+273.15;
applyBoundaryCondition(modelt,'edge',[5,6],'u',273.15); % surface temperature
applyBoundaryCondition(modelt,'edge',2,'g',0.025*ther.k); % constant moho heat flow
applyBoundaryCondition(modelt,'edge',1,'u',thermalgradient); % far field temperature
specifyCoefficients(modelt,'m',0,'d',rho*ther.cp,'c',ther.k,'a',0,'f',0);
%% 
modeln.Mesh=node2mesh(msh.Nodes,gn);
model.Mesh=node2mesh(msh.Nodes,g);
modelt.Mesh=node2mesh(msh.Nodes,gt);
%% Coeffecient
xq1 = linspace(0,0,mod.Y-1+1); % 1m grid
yq1 = linspace(0,mod.Y,mod.Y+1);
% rockstrength=1e6; % pa
stresslimitplot=((rho*9.8*yq1+atm)*0.6)'+faultstrength;
stresslowerlimitplot=((rho*9.8*yq1+atm)*0.6)';
loops=100000;
elapsedtime=0;
faultdepth0=[];
%% variablenodesize=size(model.Mesh.Nodes,2);
% initial condition
nodesize=size(model.Mesh.Nodes,2);
stressyx=zeros(nodesize,1);
stressyz=zeros(nodesize,1);
resultstressyx=createPDEResults(model,stressyx);
%% specifyCoefficients
specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',@fcoeffunction);
%% Pressure and Temperature (1-D linear field)
Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
Pfield=model.Mesh.Nodes(2,:)*rho*9.8+atm;%[Pa]
Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
%% calculate faultstress and overupperlimit
faultnodes=find(model.Mesh.Nodes(1,:)==0); % find nodes on the fault
faultnodesdepth=model.Mesh.Nodes(2,faultnodes);
faultstress=stressyx(faultnodes);
faultslip=faultstress*0;
equ=zeros(length(faultslip),1);
% sort data
[faultnodesdepth,faultstress]=sortdata(faultnodesdepth,faultstress);
ii=0;
jj=0;
interseismiccount=0;
intertime=0;
coseismiccount=0;
earthquakecount=1;
next=true; % start from interseismic
savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
timevsdepth=zeros(loops,2);
[p,e,t]=meshToPet(msh);
% [p,e,t]=meshToPet(mshg);
bmatrix=struct('Q',{},'G',{},'H',{},'R',{});
[bmatrix(1).Q,bmatrix(1).G,bmatrix(1).H,bmatrix(1).R]=assembleBoundary(model); 
Kinter=createGlobalKFGPU(p,t,'K',model.EquationCoefficients.CoefficientAssignments.c); %stiffness matrix
heating=shearheating(msh,Tfield,esyxv,esyzv);
femodel=createfemodel(modelt,Tfield,coefstruct,'new',heating);
for j=1:loops
    ii=ii+1;
    next0=next;
    next=isinterseismic(interseismiccount,faultslip,faultstress,faultnodesdepth);
    if next0~=next
        if next % switch from coseismic to interseismic
            earthquakecount=earthquakecount+1;
            coseismiccount=0;
            modeln.Mesh=msh;
            intertime=0;
            disp(max(abs(equ)))
            % frictional heating
            modelt.Mesh=node2mesh(msh.Nodes,gt);
            Tfield=addfrictionalheating(Tfield,equ,mshg);
            equ=zeros(length(faultslip),1); % reset earthquake displacement
        else % switch from interseismic to coseismic
            interseismiccount=0;    
        end
    end
    if next
        cut=0;
        coseismiccount=0;
        interseismiccount=interseismiccount+1;
        timeforcal=now;
        %% cal. stress
        % update Fdisl
        Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
        [stressyx,stressyz,elapsedtime,v,dt]=stressaccumulationp(Fdisl,stressyx,stressyz,elapsedtime,model,bmatrix,Kinter);
        %% calculate new Tfield
        % set initial conditions with old temperature field
        heating=shearheating(msh,Tfield,esyxv,esyzv);
        femodel=createfemodel(modelt,coefstruct,'update',heating);
        T1=solvelinearTDPDE(Tfield,[0,dt],femodel);
        dT=T1(:,2)'-Tfield;
        Tfield=T1(:,2)';
        % calculate the temperature field 
        MaxTdiff=max(abs(dT))*sign(dT);
        %% 
        timeforcal = (now-timeforcal)*3600*24;
        intertime=intertime+dt;
        if interseismiccount==1
            fprintf(1,'timeforcal      elapsedtime     timestep    Time since last earthquake   Max temperature change\n') ;
        end
        fprintf(+1, '%10.2f %16.2f  %11.2f  %20.2f %10.2f \n',[timeforcal,elapsedtime/(365*3600*24),dt/(365*3600*24),intertime/(365*3600*24),MaxTdiff]) ;
        %disp(['timeforcal.:' num2str(timeforcal) '/elapsedtime:' num2str(elapsedtime/(365*3600*24)) '/dt:' num2str(dt/(365*3600*24)) '\Time pased since last eq:' num2str(intertime/(365*3600*24))])
        %% plot stress
        [faultstress,faultnodesdepth]=plotfaultinfo(model,stressyx,yq1,stresslimitplot,earthquakecount,coseismiccount,elapsedtime);
        drawnow
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=v;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period='Interseismic';
    else
        interseismiccount=0;
        coseismiccount=coseismiccount+1;
       %% solve equation and update stress
        time=now;
        %% reconstruct mesh for coseismic model
        % stressyxi for stress boundary condition
        stressyxi=scatteredInterpolant(msh.Nodes(1,:)',msh.Nodes(2,:)',stressyx,'nearest');
        % cal stressfield after earthquake. update faultdepth
        [stressyx,stressyz,u]=stressdrop(stressyx,stressyz,modeln,gn);
        if ~isempty(faultdepth0)% && ~next0
            [stressyx,stressyz]=smoothstress(p,stressyx,stressyz,faultdepth0,faultdepth);
        end
        faultdepth0=faultdepth;
        time=(now-time)*3600*24;
        %% plot fault stress
        [faultstress,faultnodesdepth,faultslip]=plotfaultinfo(model,stressyx,yq1,stresslimitplot,earthquakecount,coseismiccount,elapsedtime,u);
        equ=equ+u;
        if coseismiccount==1
            fprintf(1,'timeforcal      faultdepth     maximumslip\n') ;
        end
        fprintf(+1, '%10.2f %14.2f %10.2f \n',[time,faultdepth,max(abs(faultslip)*100)]) ;
        savefield(ii).tyx=stressyx;
        savefield(ii).tyz=stressyz;
        savefield(ii).result=u;
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period=['coseismic No.' num2str(coseismiccount)];
    end
    if ii==100
        jj=jj+1;
        disp('Saving data...')
        save([savefolder filesep num2str(jj*ii,'%06.0f') '.mat']);
        disp('Data saved')
        % clear the savefield and save count
        savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
        ii=0;
    end
end
delete(poolobj)