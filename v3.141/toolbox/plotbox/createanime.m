clear
close all
addpath(genpath('toolbox'))
% datapath='/mnt/HIROTA_DISK1/170512172438/';
[dl,n]=lsdatafile;
% datafilename=dir(datapath);
datafilenum=length(dl);
figurefolder=[dl(1).folder filesep 'figure'];
delete(gcp('nocreate'))
% parpool(20);
% warning('off','all')
%% make dir
datatype={'etaeff';'esyxv';'esyzv';'tyx';'tyz';'result'};
situationtype={'In';'co'};
savefolder=cell(size(situationtype,1)*size(datatype,1),1);
for ii=1:size(situationtype,1)
    for jj=1:size(datatype,1)
        savefolder{(ii-1)*size(datatype,1)+jj}=strcat(figurefolder,filesep,datatype(jj),filesep,situationtype(ii));
        mkdir(char(savefolder{(ii-1)*size(datatype,1)+jj,1}));
    end
end
%%
i=0;
while true

    dl=lsdatafile(n);
    datafilenum=length(dl);
    if i<datafilenum
        i=i+1;
        disp([num2str(i) filesep num2str(datafilenum)])
        if ~datafilename(i).isdir
            filenum=str2double(datafilename(i).name(1:end-4))-100;
            disp(['loading file...' dl(i).name])
            load([dl(i).folder filesep dl(i).name])
            disp('file loaded')
            viscosfield=struct('etaeff',{},'esyxv',{},'esyzv',{});
            hh = waitbar(0,'Calculating data for viscous flow...');
            for k=1:length(savefield)
                waitbar(k / length(savefield))
                [viscosfield(k).etaeff,viscosfield(k).esyxv,viscosfield(k).esyzv]=stress2visco(Fdisl,savefield(k).tyx,savefield(k).tyz);
            end
            close(hh)
            %
            %         hh = waitbar(0,'ploting...');
            if faultdepth>30000
                faultdepth=30000;
            end
            if earthquakecount~=1
                climit([0,byerlee(faultdepth,'m')*1e6;min(stressyz),max(stressyz)]);
            end
            datacell1=struct2cell(viscosfield');
            datacell2=struct2cell(savefield');
            datacell=[datacell1;datacell2];
            %         period=datacell{8};
            
            for k=1:size(datacell,1)-2
                %             parfor_progress(100);
                datatype1=datatype{k};
                for j=1:size(savefield,2)
                    %                 waitbar(j / 100)
                    time=savefield(j).Elapsedtime/(365*24*3600);
                    peroid=savefield(j).period;
                    pdeplotsavedata(datacell{k,j},model,datatype1,time,j+filenum,peroid,figurefolder,climit);
                    close all
                    parfor_progress;
                end
                %             parfor_progress(0);
            end
            %
            
        end
    elseif i>=datafilenum
        disp('waiting for new file')
        pause(600);
    end
end