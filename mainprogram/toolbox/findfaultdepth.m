function [faultdepth,faultstress]=findfaultdepth(msh,stressyx,varargin)
global modc
faultnodes=find(msh.Nodes(1,:)==0); % find nodes on the fault
faultnodesdepth=msh.Nodes(2,faultnodes);
faultstress=stressyx(faultnodes);
faultstress=sortrows([faultnodesdepth',faultstress]);
stresslimitmoho=(modc.ther{1}.rho*modc.g*modc.mod.CT+modc.atm)*modc.mu;
if length(modc.ther)==1
    overlimit=faultstress(:,2)-(modc.ther{1}.rho*modc.g*faultstress(:,1)*modc.mu);
elseif length(modc.ther)==2
    overlimit=zeros(length(faultstress),1);
    overlimit(faultstress(:,1)<modc.mod.CT)=faultstress(faultstress(:,1)<modc.mod.CT,2)-(modc.ther{1}.rho*modc.g*faultstress(faultstress(:,1)<modc.mod.CT,1)*modc.mu);
    overlimit(faultstress(:,1)>=modc.mod.CT)=faultstress(faultstress(:,1)>=modc.mod.CT,2)-((modc.ther{2}.rho*modc.g*(faultstress(faultstress(:,1)>=modc.mod.CT,1)-modc.mod.CT)*modc.mu)+stresslimitmoho);
end
overupperlimit=overlimit-modc.fs;
if isempty(varargin)
    index1=find(overupperlimit>0,1,'last');
    if ~isempty(index1)
        if faultstress(index1,1)<modc.mod.CT
            overlimit(faultstress(:,1)>=modc.mod.CT)=0;
        end
    else
        overlimit(faultstress(:,1)>=modc.mod.CT)=0;
    end
else
    if varargin{1}==1
        overlimit(faultstress(:,1)>=modc.mod.CT)=0;
    end
end
% if length(modc.ther)==2
%     if ~modc.allowmantleearthquake
%         overlimit(plotdata(:,1)>modc.mod.CT)=-inf;
%     end
% end
index=find(overlimit>0);
faultdepth=faultstress(max(index)+1,1);
end
