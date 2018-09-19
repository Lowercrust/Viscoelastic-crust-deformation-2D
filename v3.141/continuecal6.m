clear
close all
global AB BA
addpath(genpath('toolbox'))
% reset(gpuDev)
gpustatus(false);
gpustatus(false,'calstressdrop')
dl=lsdatafile;
% start from save file
prompt = ['select a data file (1~' num2str(length(dl)) ')'];
filenum=input(prompt);
if exist([dl(filenum).folder filesep 'basicmodel' filesep 'basicmodel.mat'],'file')
    load([dl(filenum).folder filesep 'basicmodel' filesep 'basicmodel.mat'])
else
    load([dl(filenum).folder filesep dl(1).name])
end
disp(['loading data from' dl(filenum).folder filesep dl(filenum).name])
load([dl(filenum).folder filesep dl(filenum).name])
disp('data loaded');
if mean(stressyx)==0
    stressyx=double(savefield(end).tyx);
    stressyz=double(savefield(end).tyz);
end
savefolder=dl(1).folder;% make sure that file save at same folder
% savefolder=['..' filesep 'datasave' filesep 'v' num2str(ver) Simulationstarttime];
% savefolder=[filesep 'mnt' filesep 'disk1' filesep 'v' num2str(ver) Simulationstarttime];

clear dl
%     right=createright(msh,dt,savefield(100).tyx,savefield(100).tyz,esyxv,esyzv);
if considerheating
    savefield=struct('tyx',{},'tyz',{},'result',{},'Temp',{},'Elapsedtime',{},'period',{});
else
    savefield=struct('tyx',{},'tyz',{},'result',{},'Elapsedtime',{},'period',{});
end
savelist={'earthquakecount','elapsedtime','equ','faultdepth','faultdepth0',...
    'interseismiccount','intertime','next','next0','ii','j','jj','stressyx','stressyz'};
% if considerheating
%     [model,g,modeln,gn,modelt,gt]=createmodel(considerheating,D);
%     modeln.Mesh=node2mesh(modc.Nodes,gn);
%     model.Mesh=node2mesh(modc.Nodes,g);
%     modelt.Mesh=node2mesh(modc.Nodes,gt);
% else
%     [model,g,modeln,gn]=createmodel(considerheating);
%     modeln.Mesh=node2mesh(modc.Nodes,gn);
%     model.Mesh=node2mesh(modc.Nodes,g);
% end
ii=0;
% faultnodesdepth=sort(faultnodesdepth);
% faultnodesdepth=p(2,faultnodesindex);
% stresslimit=(modc.ther{1}.rho*modc.g*faultnodesdepth+modc.atm)*modc.mu+ones(1,length(faultnodesdepth)).*(modc.fs);
% stresslimitmoho=(modc.ther{1}.rho*modc.g*modc.mod.CT+modc.atm)*modc.mu+(modc.fs);
% stresslimit(faultnodesdepth>modc.mod.CT)=(modc.ther{2}.rho*modc.g*(faultnodesdepth(faultnodesdepth>modc.mod.CT)-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho;
% yq1=p(2,p(1,:)==0);
% [yq1,sortindex]=sortrows([yq1',stresslimit']);
% faultstress=sortrows([p(2,p(1,:)==0)',stressyx(p(1,:)==0)]);
% faultstress=faultstress(:,2);
% stressyx=double(savefield(90).tyx);
% stressyz=double(savefield(90).tyz);
% v=double(savefield(90).result);
% rej=90;
% next0=true;
% next=true;
if next0~=next
    coseismiccount=0;
    if interseismiccount==0
        gpustatus(false,'calstressdrop')
    end
    [stressyx,stressyz]=smoothstressrange(p,stressyx,stressyz,[0 0 modc.mod.X 1000],100,0.0001); %[x0 y0 x1 y1]
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
end
% delete(gcp('nocreate'))
% poolobj = parpool(4);
rej=j+1;
[A,B]=createAB(p,t);
AB=B*A;
[A,B]=createBA(p,t);
BA=B*A;
disp(num2str(earthquakecount));
tic
gpuDev=gpuDevice(1);
toc
disp('loop restarts')
for j=rej:loops
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
                if exist([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*100+ii,'%04.0f') '_data.mat'],'file')
                    error('file exist')
                else
                    disp(['Saving data to ' savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*100+ii,'%04.0f') '_data.mat'])
                    % save data struct and loop information
                    save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*100+ii,'%04.0f') '_data.mat'],'savefield',savelist{:},'-v7.3');
                    disp('Data saved')
                end
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
%         reset(gpuDev);
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
        if coseismiccount==1
            if faultdepth>modc.mod.CT
                eqtype=0; %mantle eq
            else
                eqtype=1;
            end
        end
        % stressyxi for stress boundary condition
        stressyxi=gpuArray([yq1(:,1),faultstress]);
        % cal stressfield after earthquake.
        [stressyx,stressyz,u]=stressdrop(stressyx,stressyz,modeln,gn,eqtype);
        %         reset(gpuDev);
        equ=equ+u;
        if coseismiccount~=1%~isempty(faultdepth0)% && ~next0
            [stressyx,stressyz]=smoothstress(modc.Nodes,stressyx,stressyz,faultdepth0,faultdepth);
        end
        faultdepth0=faultdepth;

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
        end
        jj=jj+1;
        disp(['Saving data to ' savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*100+ii,'%04.0f') '_data.mat'])
        %         save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*ii,'%04.0f') '.mat'],'-v7.3');
        % save data struct
        if exist([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*ii,'%04.0f') '_data.mat'],'file')
            error('file exist')
        else
            save([savefolder filesep num2str(earthquakecount,'%04.0f') '_' num2str(jj*ii,'%04.0f') '_data.mat'],'savefield',savelist{:},'-v7.3');
        end
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