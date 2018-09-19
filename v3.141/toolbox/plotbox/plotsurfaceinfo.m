
function [pdx,pdy]=plotsurfaceinfo(node,data,varargin)
p = inputParser;
addRequired(p,'node',@checksize);
addRequired(p,'data',@checkdata);
addOptional(p,'outputfig',true);
addOptional(p,'field',[]);
% addOptional(p,'range',[]);
parse(p,node,data,varargin{:})
% surfacenodes=find(); % find nodes on the fault
nodesdistance=node(1,node(2,:)==0);
if ismatrix(data) && ~isstruct(data)
    surfacedata=data(node(2,:)==0);
    plotdata=sortrows([nodesdistance',surfacedata]);
    pdx=plotdata(:,1);
    pdy=plotdata(:,2);
    if p.Results.outputfig
        plot(pdx,pdy);
    end
elseif isstruct(data)
    cl=[];
    names= fieldnames(data);
    for i=1:length(names)
        if strcmp(p.Results.field,names{i})
            cl=i;
            break
        end
    end
    data=struct2cell(data);
%     pdx=cell(length(data),1);
    pdy=cell(length(data),1);
    for i=1:length(data)
        datanow=data{cl,1,i};
        surfacedata=datanow(node(2,:)==0);
        plotdata=sortrows([nodesdistance',surfacedata]);
%         disp(i)
%         pdx{i}=plotdata(:,1);
        pdy{i}=plotdata(:,2);
    end
end
pdx=plotdata(:,1);
end
function result=checksize(input)
    if ismatrix(input) && size(input,1)==2 
        result=true;
    else
        result=false;
    end
end
function result=checkdata(input)
    if ismatrix(input) || isstruct(input)
        result=true;
    else
        result=false;
    end
end