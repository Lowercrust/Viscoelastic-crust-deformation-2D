clear
close all
disp('Nvidia GPGPU with at least 3GB gpu memory is required to run this program,')
pause;
addpath(genpath('toolbox'))
gpustatus(false)
gpustatus(false,'calstressdrop')
global modc faultdepth stressyxi gpunum Tfield0 fi%sliplength
gpunum=gpuDeviceCount;
formatOut = 'yymmddHHMMSS';
Simulationstarttime=datestr(now,formatOut);
ver=3.141;
addpath(genpath('toolbox'))
savefolder=[filesep 'mnt' filesep 'disk1' filesep 'v' num2str(ver) Simulationstarttime];
% savefolder=['..' filesep 'datasave' filesep 'v' num2str(ver) Simulationstarttime];
modc=constants;
prompt = 'Consider heating ?(y or n)';
considerheating=input(prompt,'s');
% initial faultdepth, just for creating model geometry.
faultnodesindex=find(modc.Nodes(1,:)==0);
faultnodesy=modc.Nodes(2,faultnodesindex);
faultdepth=min(faultnodesy(faultnodesy>0));
Tb=modc.mod.Tb;
Tfun=@(D,Tb) geothermal(D)-Tb;
Tfun1=@(D) Tfun(D,Tb);
D=fsolve(Tfun1,1000);
[Tfield,A0]=geothermal(modc.Nodes(2,:),D);%[K]
mkdir(savefolder)
% yes for consider heating
if strcmp(considerheating,'y')
    considerheating=true;
    disp('creating model')
    [model,g,modeln,gn,modelt,gt]=createmodel(considerheating,D);
    disp('preparing femodel for heat flow equation...')
    femodel=createfemodel(modelt,'new');
    save([savefolder filesep 'femodel.mat'],'femodel','-v7.3');
    savefield=struct('tyx',{},'tyz',{},'result',{},'Temp',{},'Elapsedtime',{},'period',{});
else
    considerheating=false;
    [model,g,modeln,gn]=createmodel(considerheating);
    savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
end
[p,e,t]=meshToPet(modeln.Mesh);

disp('preparing femodel for stress equation (interseismic)...')
femodelinter=createfemodel(model,'new');
% save femodel
save([savefolder filesep 'femodelinter.mat'],'femodelinter','-v7.3');
%% Pressure and Temperature (1-D linear field)

% Tfield=Tfield;
if size(modc.rhe,1)==2
    Pfield=zeros(1,length(modc.Nodes));
    fi=cell(2,1);
    fi{1}=find(modc.Nodes(2,:)<=modc.mod.CT);
    fi{2}=find(modc.Nodes(2,:)>modc.mod.CT);
    Pfield(1,fi{1})=modc.Nodes(2,fi{1})*modc.ther{1}.rho*modc.g+modc.atm;%[Pa]
    Pfield(1,fi{2})=modc.Nodes(2,fi{2})*modc.ther{2}.rho*modc.g+modc.atm;%[Pa]
    DH=zeros(1,length(modc.Nodes));
    % Decay heating
    DH(fi{1})=(A0*exp(-modc.Nodes(2,fi{1})/D));
    DH=pdeintrpgpu(p,t,DH');
else
    fi=cell(1,1);
    fi{1}=(1:length(modc.Nodes))';
    Pfield=modc.Nodes(2,:)*modc.ther.rho*modc.g+modc.atm;%[Pa]
end
Fdisl=FPT(Pfield,Tfield,modeln.Mesh,'disl')'; % [MPa^n*s]
Tfield0=Tfield;%save for initial temperature field
%% out put basic information for the data set
createmodelinfo(Simulationstarttime)
%% create femodel for interseismic model
%% Coeffecient
loops=100000;
elapsedtime=0;
faultdepth0=[];
%% variablenodesize=size(model.Mesh.Nodes,2);
% initial condition
nodesize=size(modc.Nodes,2);
stressyx=zeros(nodesize,1);
stressyz=zeros(nodesize,1);
%% calculate faultstress and overupperlimit
faultnodes=find(modc.Nodes(1,:)==0); % find nodes on the fault
faultnodesdepth=modc.Nodes(2,faultnodes);
faultstress=stressyx(faultnodes);
faultslip=faultstress*0;
% sort data
% [faultnodesdepth,faultstress]=sortdata(faultnodesdepth,faultstress);
ii=0;
jj=0;
interseismiccount=0;
intertime=0;
coseismiccount=0;
earthquakecount=1;
next=true; % start from interseismic
calinfo=struct('count',{},'timeforcal',{},'et',{},'dt',{},'intertime',{},'dT',{},'dT0',{},'faultdepth',{},'maxslip',{});
calinfo(1).intertime=0;
equ=zeros(length(modc.Nodes),1);
%
faultnodesdepth=p(2,faultnodesindex);
stresslimit=(modc.ther{1}.rho*modc.g*faultnodesdepth+modc.atm)*modc.mu+ones(1,length(faultnodesdepth)).*(modc.fs);
stresslimitmoho=(modc.ther{1}.rho*modc.g*modc.mod.CT+modc.atm)*modc.mu+(modc.fs);
stresslimit(faultnodesdepth>modc.mod.CT)=(modc.ther{2}.rho*modc.g*(faultnodesdepth(faultnodesdepth>modc.mod.CT)-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho;
yq1=p(2,p(1,:)==0);
[yq1,sortindex]=sortrows([yq1',stresslimit']);
faultstress=sortrows([p(2,p(1,:)==0)',stressyx(p(1,:)==0)]);
faultstress=faultstress(:,2);
gpuDev=gpuDevice;
eqtype=[];
%%
% save basicmodel
save([savefolder filesep 'basicmodel.mat']);
% list of non constants that change with loop
savelist={'earthquakecount','elapsedtime','equ','faultdepth','faultdepth0','interseismiccount','intertime','next','next0','ii','j','jj','stressyx','stressyz'};
disp('loop starts')
% ii: step counter(1,100) jj:file counter 
for j=1:loops
    ii=ii+1; 
    next0=next; 
    next=isinterseismic(interseismiccount,faultslip,faultstress,yq1(:,1)',yq1(:,2)');
    if next0~=next
        if next % switch from coseismic to interseismic
            coseismiccount=0;
            if interseismiccount==0
                gpustatus(false,'calstressdrop')
            end
            calinfo=struct('count',{},'timeforcal',{},'et',{},'dt',{},'intertime',{},'dT',{},'dT0',{},'faultdepth',{},'maxslip',{});
            calinfo(1).intertime=0;
            % smooth stress near surface after each earthquake
            [stressyx,stressyz]=smoothstressrange(p,stressyx,stressyz,[0 0 modc.mod.X 1000],100,0.0001); %[x0 y0 x1 y1]
            if considerheating
                Tfield=addfrictionalheating(Tfield,equ,modelt.Mesh);
            end
            % save data after each earthquake
            if ~exist(savefolder,'dir')
                mkdir(savefolder)
            end
            if ii~=1 % if ii=1, data already saved
                disp('Saving data...')
                % save data struct and loop information
                save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*100+ii,'%04.0f') '_data.mat'],'savefield',savelist{:},'-v7.3');
                disp('Data saved')
            end
            % reset
            equ=zeros(length(modc.Nodes),1);
            earthquakecount=earthquakecount+1;
            ii=1;
            jj=0;
            % reset savefield
            if considerheating
                savefield=struct('tyx',{},'tyz',{},'result',{},'Temp',{},'Elapsedtime',{},'period',{});
            else
                savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
            end
        else % switch from interseismic to coseismic
            if coseismiccount==0
                gpustatus(true,'calstressdrop')
            end
            interseismiccount=0;
            calinfo=struct('count',{},'timeforcal',{},'et',{},'dt',{},'intertime',{},'dT',{},'dT0',{},'faultdepth',{},'maxslip',{});
            calinfo(1).dt=0;
            eqtype=[];
        end
    end
    %% iteration
    if next
        time=now;
        cut=0;
        coseismiccount=0;
        interseismiccount=interseismiccount+1;
        %% calculate new stress field
        [stressyx,stressyz,etaeff,esyxv,esyzv,elapsedtime,v,dt]=stressaccumulationp(Fdisl,stressyx,stressyz,elapsedtime,model,femodelinter);
%         disp(gpuDev.FreeMemory)
        reset(gpuDev);
         %% calculate new Tfield
        % set initial conditions with old temperature field
        if considerheating
            heating=shearheating(modelt.Mesh,Pfield,Tfield,esyxv,esyzv)+DH;
            femodel=createfemodel(modelt,'update',heating,femodel);
            [Tfield,MaxdiffdT,MaxdiffdT0]=solvelinearTDPDE(Tfield,[0,dt],femodel);
            % calculate the temperature field
            Fdisl=FPT(Pfield,Tfield,modelt.Mesh,'disl')'; % [MPa^n*s] update Fdisl
            savefield(ii).Temp=single(Tfield);
            calinfo(1).dT=MaxdiffdT;
            calinfo(1).dT0=MaxdiffdT0;
        end
        %% updata and display status of calculation 
        calinfo(1).timeforcal = (now-time)*3600*24;
        calinfo(1).intertime=calinfo.intertime+dt;
        calinfo(1).elapsedtime=elapsedtime;
        calinfo(1).dt=dt;
        calinfo(1).count=interseismiccount;
        dispcalculationinfo(calinfo)
        %% plot stress
        [faultstress,faultnodesdepth]=plotfaultinfo(model,stressyx,yq1(:,1),yq1(:,2),earthquakecount,coseismiccount,elapsedtime);
        savefield(ii).tyx=single(stressyx);
        savefield(ii).tyz=single(stressyz);
        savefield(ii).result=single(v);
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period='Interseismic';
    else
        
        interseismiccount=0;
        coseismiccount=coseismiccount+1;
       %% solve equation and update stress
        time=now;
        % stressyxi for stress boundary condition
        stressyxi=gpuArray([yq1(:,1),faultstress]);
        % cal stressfield after earthquake.
        [stressyx,stressyz,u]=stressdrop(stressyx,stressyz,modeln,gn,eqtype);
        reset(gpuDev);
        equ=equ+u;
        if coseismiccount~=1%~isempty(faultdepth0)% && ~next0
            [stressyx,stressyz]=smoothstress(modc.Nodes,stressyx,stressyz,faultdepth0,faultdepth);
        end
        faultdepth0=faultdepth;
        if coseismiccount==1
            if faultdepth>modc.mod.CT
                eqtype=0; %mantle eq
            else
                eqtype=1;
            end
        end
        time=(now-time)*3600*24;
        calinfo(1).timeforcal = time;
        calinfo(1).faultdepth=faultdepth;
        calinfo(1).maxslip=max(abs(equ));
        calinfo(1).count=coseismiccount;
        dispcalculationinfo(calinfo)
        %% plot fault stress
        [faultstress,faultnodesdepth,faultslip]=plotfaultinfo(model,stressyx,yq1(:,1),yq1(:,2),earthquakecount,coseismiccount,elapsedtime,u);
        drawnow
        savefield(ii).tyx=single(stressyx);
        savefield(ii).tyz=single(stressyz);
        savefield(ii).result=single(u);
        savefield(ii).Elapsedtime=elapsedtime;
        savefield(ii).period=['coseismic No.' num2str(coseismiccount)];
    end
    if ii==100
        if ~exist(savefolder,'dir')
            mkdir(savefolder)
            disp([savefolder ' has been created']);
        end
        jj=jj+1;
        disp(['Saving data to ' savefolder])
        %         save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*ii,'%04.0f') '.mat'],'-v7.3');
        % save data struct
        save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*ii,'%04.0f') '_data.mat'],'savefield',savelist{:},'-v7.3');
        disp('Data saved')

        % clear the savefield and save count
        if considerheating
            savefield=struct('tyx',{},'tyz',{},'result',{},'Temp',{},'Elapsedtime',{},'period',{});
        else
            savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
        end
        ii=0;
    end
end
delete(poolobj)