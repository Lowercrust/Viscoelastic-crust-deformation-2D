% out put mat file show gpu status
% true or false for request;
% true: gpu is busy
function gpustatus(request,varargin)
% gpustat=request;
cname=getComputerName;
if ~strcmp(cname,'red.seis.net')
if ~isempty(varargin)
    cname=[cname varargin{1}];
end
cname=[cname '.gpustat'];
if request
    while exist(cname,'file')==2
%         pause(0.5)
    end
    save(cname,'request')
else
    while exist(cname,'file')==2
        delete(cname)
    end
end
end
end
