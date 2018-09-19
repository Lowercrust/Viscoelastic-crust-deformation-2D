% Create rectangular which has n node point on boundary
% fault has been divided into n-3 sections
% n is larger than or equal to 5
function [x,y]=rect(bs,s)
global mod
X=mod.X;
Y=mod.Y;
nbs=4; % number of boundary segments
% startpointlist=[0 slipregion(cutsidx+1)];
% endpointlist=[slipregion(cutsidx-1) slipregion(end)];
if nargin==0,
    x=nbs; % number of boundary segments
    return
end
d=diag([0,1,0,1])*ones(4,nbs);
d(end-1:end,end-1:end)=[1,1;0,0];
% d=[
%     0 0 0 0 0 0 % start parameter value
%     1 1 1 1 1 1 % end parameter value
%     0 0 0 0 1 1 % left hand region
%     1 1 1 1 0 0 % right hand region
%     ];
bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
    error(message('pde:lshapeg:InvalidBs'))
end

if nargin==1,
    x=d(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
    error(message('pde:lshapeg:SizeBs'));
end
if ~isempty(s),
    % boundary segment 1 to nbs-3
    
    ii=find(bs==1);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[0 Y],s(ii));
    end
    
    % boundary segment nbs-2
    ii=find(bs==2);
    if ~isempty(ii)
        x(ii)=interp1([d(1,nbs-2),d(2,nbs-2)],[0 X],s(ii));
        y(ii)=interp1([d(1,nbs-2),d(2,nbs-2)],[Y Y],s(ii));
    end
    
    % boundary segment nbs-1
    ii=find(bs==3);
    if ~isempty(ii)
        x(ii)=interp1([d(1,nbs-1),d(2,nbs-1)],[X X],s(ii));
        y(ii)=interp1([d(1,nbs-1),d(2,nbs-1)],[0 Y],s(ii));
    end
    
    % boundary segment nbs
    ii=find(bs==4);
    if ~isempty(ii)
        x(ii)=interp1([d(1,nbs),d(2,nbs)],[0 X],s(ii));
        y(ii)=interp1([d(1,nbs),d(2,nbs)],[0 0],s(ii));
    end
end