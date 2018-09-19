% function [esyxvcell,esyzvcell,t]=viscousrate
clear
close all
% global modc fi
dl=lsdatafile;
load([dl(1).folder filesep 'basicmodel' filesep 'basicmodel']);
load(['matfile' filesep 'stress0.mat']);
[etaeff0,esyxv,esyzv]=stress2visco(Fdisl,model.Mesh,stressyx0,stressyz0);
% [pindex,x1,y1]=grid2pindex(x,y,p);
t=cell(length(dl),1);
savefolder=[dl(1).folder filesep 'viscousrate'];
mkdir(savefolder);
for i=41:length(dl)
    disp(['loading data' dl(1).folder filesep dl(i).name]);
    load([dl(1).folder filesep dl(i).name]);
    esyxvcell=cell(100,1);
    time=zeros(length(savefield),1);
    if strcmp(dl(i).name(11:end),'lineardata.mat')
        for j=1:length(savefield)
            if strcmp(savefield(j).period,'Interseismic')
                [esyxv,esyzv]=stress2viscolinear(stressyx0,stressyz0,etaeff0);
                esyxvcell{j}=esyxv;
                %             esyzvcell{j}=esyzv;
                time(j)=savefield(j).Elapsedtime;
            end
        end
    else
        %     esyzvcell=cell(100,1);
        for j=1:length(savefield)
            if strcmp(savefield(j).period,'Interseismic')
                [etaeff,esyxv,esyzv]=stress2visco(Fdisl,model.Mesh,savefield(j).tyx,savefield(j).tyz);
                esyxvcell{j}=esyxv;
                %             esyzvcell{j}=esyzv;
                time(j)=savefield(j).Elapsedtime;
            end
        end
    end
    t{i}=time;
    save([savefolder filesep num2str(i,'%04.0f') ...
        '_viscousrate.mat'],'esyxvcell','t','-v7.3');
end