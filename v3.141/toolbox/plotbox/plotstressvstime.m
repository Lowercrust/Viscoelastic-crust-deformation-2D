% plot stress changes at a specific location
% multiple locations can be specified [x1,x2,x3],[y1,y2,y3]
function [datacell,timedata]=plotstressvstime(x,y)
dl=lsdatafile;
filenum=length(dl);
load([dl(end).folder filesep dl(1).name]);
% find nearest node
ni=zeros(length(x),1);
data=zeros(length(savefield),filenum);
datacell=cell(length(x),1);
timedata=data;
for i=1:length(x)
    distance=hypot(modc.Nodes(1,:)-x(i),modc.Nodes(2,:)-y(i));
    [~,ni(i)]=min(distance);
    disp(modc.Nodes(:,ni(i)));
    datacell{i}=zeros(length(savefield),filenum);
end
h = waitbar(0,'Please wait...');
for i=1:filenum
    waitbar(i/filenum,h,sprintf('%3.0f',i))
    load([dl(end).folder filesep dl(i).name],'savefield');
    if length(savefield)>1
        savefieldtable=struct2table(savefield);
        timedata(1:length(savefield),i)=table2array(savefieldtable(:,4));
        for k=1:length(ni)
            for j=1:length(savefield)
                datacell{k}(j,i)=savefield(j).tyx(ni(k));
            end
        end
    else
        timedata(1,i)=savefield.Elapsedtime;
        for k=1:length(ni)
            datacell{k}(1,i)=savefield.tyx(ni(k));
        end
    end
end
end