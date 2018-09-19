%interpolate data from mshold to mshnew
function varargout=meshchange(mshold,mshnew,isoutputfaultinfo,varargin)
oldx=mshold.Nodes(1,:)';
oldy=mshold.Nodes(2,:)';
newx=mshnew.Nodes(1,:)';
newy=mshnew.Nodes(2,:)';

varargout=varargin;
datacell=varargin;
parfor i=1:length(varargin)
    %disp(i)
%     disp(size(varargin{i}))
    datacell{i}=scatteredInterpolant(oldx,oldy,varargin{i});
    varargout{1,i}=datacell{i}(newx,newy);
end
% interpolate data from interseismic to coseismic
if isoutputfaultinfo
    faultnodes=find(mshnew.Nodes(1,:)==0); % find nodes on the fault
    varargout{length(varargin)+1}=mshnew.Nodes(2,faultnodes); % faultnodesdepth
    varargout{length(varargin)+2}=stressyx(faultnodes);% faultstress
end
end