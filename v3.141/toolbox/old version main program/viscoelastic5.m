clear
close all
global mod rhe yieldstress ther
global R rho right result resultn pointlist overlimit region1 resultstressyx rockstrength %sliplength
formatOut = 'yymmddHHMMSS';
Simulationstarttime=datestr(now,formatOut);
savefolder=['H:\datasave\' Simulationstarttime];
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
% Heatflow model
modelh=model;
%% Boundary Conditions
% Boundary Conditions for mechanical equation 
applyBoundaryCondition(model,'edge',2,'u',v0);
applyBoundaryCondition(model,'edge',4,'u',0);
% Boundary Conditions for heat flow equation
applyBoundaryCondition(modelh,'edge',1,'u',273.15);
applyBoundaryCondition(modelh,'edge',3,'g',0.025);
%% Coeffecients
faultstress0=zeros(mod.ny,1);
xq1 = linspace(0,0,mod.Y+1); % 1m grid
yq1 = linspace(0,mod.Y,mod.Y+1);
maxoverlimit=0;
rockstrength=1e6;
stresslimit=(rho*9.8*yq1*0.6)'+rockstrength;
stresslowerlimit=(rho*9.8*yq1*0.6)';
earthquakes=1000;
Events=struct;
poolobj = parpool(16);
overupperlimit=0;
faultdepth=0;
dt=3600*24*365*1000;%[s]
dt0=dt;
elapsedtime=0;
faultstress=zeros(size(yq1,2),1);
faultslip=faultstress;
for j=1:earthquakes
    Earthquake=struct('InterseismicStress',{},'Elapsedtime',{},'CoseismicStress',{},'CoseismicSlip',{});
    %(1) stressyx (2)stressyz (3)total shear strain rate yx (4) total shear
    %strain rate yz (5) interseismic or coseismic
    savefield=struct('tyx',{},'tyz',{},'result',{},'period',{});
    if exist('interseismicmodel','var')
        model=interseismicmodel;
    end
    generateMesh(model,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
    generateMesh(modelh,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
    model.Mesh=pointrefine(g,model.Mesh,0,faultdepth,5000,1);
    modelh.Mesh=model.Mesh;
    %Coefficients for heat flow equation, f is shear heating
    specifyCoefficients(modelh,'m',0,'d',rho*ther.cp,'c',ther.k,'a',0,'f',0);
    

    %% variables
    if j==1
        %% Pressure and Temperature (1-D linear field)
        Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
        Pfield=model.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
        Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
        setInitialConditions(modelh,Tfield);
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
    end
    i=0;
    disp(['earthquake No.' num2str(j)])
    %% Interseismic faultstress<stresslimit
    while isempty(overupperlimit(overupperlimit>0))
        i=i+1;
        if exist('modeln','var') && i==1
            % interpolate stressfield after earthquake to interseismic
            % model. modeln¨model.
            stressyx=meshchange(model,modeln,stressyx);
            stressyz=meshchange(model,modeln,stressyz);
            if ~isempty(find(isnan(stressyx),1)) || ~isempty(find(isnan(stressyz),1))
                pause
            end
            stress=sqrt(stressyx.^2+stressyz.^2);
        end
        %% stress ¨@effective viscosity and viscous shear strain rate
        [etaeff,esyxv,esyzv]=stress2visco(stress,Fdisl,stressyx,stressyz);
        %disp(['size of Meshnodes ' num2str(length(model.Mesh.Nodes))])
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
        specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',@fcoeffunction);
        righttest=-interpolateSolution(right,region1.x,region1.y);
        %     disp(mean(righttest)/mean(f1))
        %% solve the equation for total shear strain rate
        %disp('Solving velocity field')
        tic
        result = solvepde(model);
        timetosolve=toc;
        v = result.NodalSolution; %[m s^-1]
        % Total shear strain rate
        esyx=result.XGradients;% [s^-1]
        esyz=result.YGradients;% [s^-1]
        %% stress accumulation from t_{0} to t_{0}+ƒ¢t
        stressyx=stressyx+(esyx-esyxv)*rhe.G*dt; %[Pa] stress of previous time step; elas stress accumulation; stress drop due to earthquake; -stressdropyx
        stressyz=stressyz+(esyz-esyzv)*rhe.G*dt; %[Pa]
        stress=sqrt(stressyx.^2+stressyz.^2); %[Pa]
        resultstressyx=createPDEResults(model,stressyx);
        elapsedtime=elapsedtime+dt;
        
        resultstressyx1=resultstressyx;
        tic
        parfor ii=1:size(yq1,2)
            faultstress(ii) = interpolateSolution(resultstressyx1,xq1(ii),yq1(ii));
        end
        timetointerploate=toc;
        disp(['Totel Elapsedtime ' num2str(elapsedtime/(365*24*3600)) '/Time to solve ' num2str(timetosolve) '/Time to interploate ' num2str(timetointerploate)])
        %% save
        Earthquake(i).InterseismicStress=faultstress;
        Earthquake(i).Elapsedtime=elapsedtime/(365*24*3600);
        savefield(i).tyx=stressyx;
        savefield(i).tyz=stressyz;
        savefield(i).result=result;
        savefield(i).period='Interseismic';
        %% plot
        plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k')
        if max(faultstress)>0
            ylim([0,max(faultstress)*2]);
            xlim([0,max(faultstress)*2/(rho*9.8)]);
        end
        title(['Interseismic' num2str(j) ',' num2str(i)])
        drawnow
        % cal. overupperlimit
        overupperlimit=faultstress-stresslimit;
    end
    interseismicmodel=model;
    %% Eartbquakes faultstress>stresslimit
    jj=0;
    u=[];
    totalfaultslip=zeros(size(yq1,2),1);
    while ~isempty(overupperlimit(overupperlimit>0))
        jj=jj+1;
        disp(['Event No.' num2str(jj)]);
        % create model
        gn=@npointrect;
        % find fault depth
        overlimit=faultstress-stresslowerlimit;
        [pointlist,nbs]=findsection(overlimit);
        faultdepth=pointlist(end-1);
        modeln = createpde;
        geometryFromEdges(modeln,gn);
        % apply boundary condition
        if nbs==4
            applyBoundaryCondition(modeln,'edge',3,'u',0);
            applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','off');
        else
            if rem(nbs,2)==1
                applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
                if nbs-4>2
                    applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'u',0);
                end
                applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'g',@bcneumann,'Vectorized','off');
            elseif rem(nbs,2)==0
                applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
                if nbs-4>1
                    applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'u',0);
                end
                applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'g',@bcneumann,'Vectorized','off');
            end
        end
        % specifyCoefficients
        specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
        u0=u;
        % ‡@ solve equation (adaptmesh)
        [u,mesh,exitflag]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',30000,'ngen',inf); % coseismic slip
        if exitflag==1
            modeln=model;
            u=u0;
            break
        end
        resultn=createPDEResults(modeln,u);
        stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake yx
        stressdropyz=rhe.G*resultn.YGradients; % yz
        stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2); % maximum shear stress drop
        % ‡A interpolate old stress field to new meshgrid (adaptmesh)
        tic
        stressyx=meshchange(modeln,model,stressyx);
        toc
        stressyz=meshchange(modeln,model,stressyz);
        toc
        % new model become old model
        model=modeln;
        % ‡B update stressfield ¨ new stress field
        stressyx=stressyx-stressdropyx;
        stressyz=stressyz-stressdropyz;
        stress=sqrt(stressyx.^2+stressyz.^2);
        % plot fault stress
        resultstressyx=createPDEResults(modeln,stressyx);
        resultstressyx1=resultstressyx;
        resultn1=resultn;
        tic

        parfor ii=1:length(yq1)
            faultstress(ii) = interpolateSolution(resultstressyx1,xq1(ii),yq1(ii));
            faultslip(ii)=interpolateSolution(resultn1,xq1(ii),yq1(ii));
        end
        totalfaultslip=totalfaultslip+faultslip;
        toc
        stresslowerlimit=(rho*9.8*yq1*0.6)';
        stresslimit=(rho*9.8*yq1*0.6+1e6)';
        %% save
        Earthquake(jj).CoseismicStress=faultstress;
        Earthquake(jj).CoseismicSlip=faultslip;
        savefield(i).tyx=stressyx;
        savefield(i).tyz=stressyz;
        savefield(i).result=resultn;
        savefield(i).period='Coseismic, stressfield after each events';
        %% plot
        plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k',yq1,-faultslip*1e9,'g')
%         if max(faultstress)>0
            ylim([0,max(faultstress)*2]);
            xlim([0,max(faultstress)*2/(rho*9.8)]);
%         end
        title(['Earthquake No. ' num2str(j) 'Event' num2str(jj)])
        drawnow
        overupperlimit=faultstress-stresslimit;
    end
    %% initial temperature field for next interseismic period
    
    %%
    EqNo=['Earthquake' num2str(j)];
    Events.(EqNo)=Earthquake;
    save([savefolder '\' EqNo '.mat'],'savefield');
end
delete(poolobj)