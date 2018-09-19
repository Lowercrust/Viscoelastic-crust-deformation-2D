clear
close all
global mod rhe yieldstress
global R rho right result resultn pointlist overlimit region1 resultstressyx rockstrength%sliplength
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
%% Boundary Conditions
applyBoundaryCondition(model,'edge',2,'u',v0);
applyBoundaryCondition(model,'edge',4,'u',0);
%% Coeffecient
faultstress0=zeros(mod.ny,1);
% xq = linspace(0,0,mod.ny); % 50m grid
% yq = linspace(0,mod.Y,mod.ny);
xq1 = linspace(0,0,mod.Y+1); % 1m grid
yq1 = linspace(0,mod.Y,mod.Y+1);
maxoverlimit=0;
rockstrength=1e6;
stresslimit=(rho*9.8*yq1*0.6)'+rockstrength;
% stresslimit=(byerlee(yq,'m')*1e6)';
stresslowerlimit=(rho*9.8*yq1*0.6)';
earthquakes=1000;
Events=struct;
% Events=cell(earthquakes,1); %recording the time and depth of each events
% savefaultstress=cell(10000,3);
% savenum=0;
poolobj = parpool(16);
overupperlimit=0;
% h1=msgbox(num2str(elapsedtime/(365*24*3600)),'elapsedtime');
faultdepth=0;
dt=3600*24*365*1000;%[s]
dt0=dt;
elapsedtime=0;
for j=1:earthquakes
    Earthquake=struct('InterseismicStress',{},'Elapsedtime',{},'CoseismicStress',{},'CoseismicSlip',{});
    %     for i=1:loops
    if exist('interseismicmodel','var')
        model=interseismicmodel;
    end
    generateMesh(model,'Hmax',500,'GeometricOrder','linear','Jiggle','on','MesherVersion','R2013a');
    %model.Mesh=regionalrefinemesh(g,model.Mesh,0,0,1000,1000,1);
    model.Mesh=pointrefine(g,model.Mesh,0,faultdepth,5000,1);
    %
    c=1;
    a=0;
    %         right=createPDEResults(model,zeros(nodesize,1));
    %% Pressure and Temperature
    Tfield=model.Mesh.Nodes(2,:)*mod.dT+273.15;%[K]
    Pfield=model.Mesh.Nodes(2,:)*rho*9.8;%[Pa]
    Fdisl=FPT(Pfield,Tfield,'disl')'; % [MPa^n*s]
    %% preset parameters
    if j==1
        nodesize=size(model.Mesh.Nodes,2);
        stressyx0=zeros(nodesize,1);
        stressyz0=zeros(nodesize,1);
        stress0=zeros(nodesize,1);
        stressyx=zeros(nodesize,1);
        stressyz=zeros(nodesize,1);
        stress=zeros(nodesize,1);
        esyxv0=zeros(nodesize,1);
        esyzv0=zeros(nodesize,1);
        esyxv=zeros(nodesize,1);
        esyzv=zeros(nodesize,1);
        stressdropyx=zeros(nodesize,1);
        stressdropyz=zeros(nodesize,1);
    end
    i=0;
    disp(['earthquake No.' num2str(j)])
    while isempty(overupperlimit(overupperlimit>0))
        i=i+1;
        %% stress Å®Å@effective viscosity and viscous shear strain rate
        if exist('modeln','var') && i==1
            stressyx=meshchange(model,modeln,stressyx);
            stressyz=meshchange(model,modeln,stressyz);
            if ~isempty(find(isnan(stressyx),1)) || ~isempty(find(isnan(stressyz),1))
                pause
            end
            stress=sqrt(stressyx.^2+stressyz.^2);
        end
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
        specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',@fcoeffunction);
        righttest=-interpolateSolution(right,region1.x,region1.y);
        %     disp(mean(righttest)/mean(f1))
        %% solve the equation for total shear strain rate
        %disp('Solving velocity field')
        tic
        result = solvepde(model);
        timetosolve=toc;
        v = result.NodalSolution; %[m s^-1]
        %         [v,mesh,code]=generateadaptMesh(model,'MesherVersion','R2013a','maxt',50000,'ngen',inf);
        %         result=createPDEResults(model,v);
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
        resultstressyx=createPDEResults(model,stressyx);
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
        disp(['Totel Elapsedtime ' num2str(elapsedtime/(365*24*3600)) '/Time to solve ' num2str(timetosolve) ])
        resultstressyx1=createPDEResults(model,stressyx);
        parfor ii=1:length(yq1)
            faultstress(ii) = interpolateSolution(resultstressyx1,xq1(ii),yq1(ii));
        end
        %         if ~isempty(find(faultstress<0,1))
        %             pause
        %         end
        Earthquake(i).InterseismicStress=faultstress;
        Earthquake(i).Elapsedtime=elapsedtime/(365*24*3600);
        %         savenum=savenum+1;
        %         savefaultstress{savenum,1}=faultstress';
        %         savefaultstress{savenum,2}=0;
        plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k')
        if max(faultstress)>0
            ylim([0,max(faultstress)*2]);
            xlim([0,max(faultstress)*2/(rho*9.8)]);
        end
        title(['Interseismic' num2str(j) ',' num2str(i)])
        drawnow
        %         anime1(i+k) = getframe(gcf);
        %         stressyx0=stressyx;
        %         stressyz0=stressyz;
        %         stress0=stress;
        %         if max(stress>byerlee(15000,'m')*1e6)
        overupperlimit=faultstress'-stresslimit;
        %         if ~isempty(overupperlimit(overupperlimit>0))
        %             k=k+i;
        %             %             break
        %         end
    end
    interseismicmodel=model;
    %% Eartbquakes
    %     for jj=1:1000
    jj=0;
    %     overupperlimit
    u=[];
    while ~isempty(overupperlimit(overupperlimit>0))
        jj=jj+1;
        disp(['Event No.' num2str(jj)]);
        % create model
        gn=@npointrect;
        % find fault depth
        overlimit=faultstress'-stresslowerlimit;
        % if isempty(pointlist)
        [pointlist,nbs]=findsection(overlimit);
        faultdepth=pointlist(end-1);
        % else
        %     pointlist0=pointlist(2:end-1);
        %     [pointlist,~]=findsection(overlimit);
        %     pointlist=sort([pointlist0 pointlist]);
        %     nbs=length(pointlist)+2;
        %     disp(pointlist)
        % end
        % pointlist=[0 35000];
        % nbs=4;
        modeln = createpde;
        geometryFromEdges(modeln,gn);
        % pdegplot(modeln, 'edgeLabels', 'on');
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
        specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
        % [p1,e1,t1]=meshToPet(model.Mesh);
        % calculate displacement and stressdrop
        %         lengthu=length(u);
        u0=u;
        % solve equation
        [u,mesh,code]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',30000,'ngen',inf); % coseismic slip
        if code==1
            modeln=model;
            u=u0;
            break
        end
        resultn=createPDEResults(modeln,u);
        stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake yx
        stressdropyz=rhe.G*resultn.YGradients; % yz
        stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2); % maximum shear stress
        % interpolate stress to adaptmesh
        tic
        stressyx=meshchange(modeln,model,stressyx);
        toc
        stressyz=meshchange(modeln,model,stressyz);
        toc
        model=modeln;
        % update stressfield
        stressyx=stressyx-stressdropyx;
        stressyz=stressyz-stressdropyz;
        stress=sqrt(stressyx.^2+stressyz.^2);
        % plot fault stress
        resultstressyx=createPDEResults(modeln,stressyx);
        resultstressyx1=resultstressyx;
        parfor ii=1:length(yq1)
            faultstress(ii) = interpolateSolution(resultstressyx1,xq1(ii),yq1(ii));
        end
        stresslowerlimit=(rho*9.8*yq1*0.6)';
        stresslimit=(rho*9.8*yq1*0.6+1e6)';
        faultslip=interpolateSolution(resultn,xq1,yq1);
        Earthquake(jj).CoseismicStress=faultstress;
        Earthquake(jj).CoseismicSlip=faultslip;
        plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k',yq1,-faultslip*1e9,'g')
        if max(faultstress)>0
            ylim([0,max(faultstress)*2]);
            xlim([0,max(faultstress)*2/(rho*9.8)]);
        end
        title(['Earthquake No. ' num2str(j) 'Event' num2str(jj)])
        drawnow
        overupperlimit=faultstress'-stresslimit;
        
        %         if isempty(overupperlimit(overupperlimit>0))
        %             break
        %         end
    end
    EqNo=['Earthquake' num2str(j)];
    Events.(EqNo)=Earthquake;
end
delete(poolobj)