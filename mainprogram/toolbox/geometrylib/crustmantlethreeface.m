function [x,y]=crustmantlethreeface(bs,s)
global modc faultdepth 
if faultdepth==0
    depthextend=100;
else
    depthextend=0;
end
if isempty(modc.fw)
    width=0.5;
else
    width=modc.fw;
end
X=modc.mod.X;
b=modc.mod.CT;
Y=modc.mod.Y;
nbs=11;
if nargin==0
    x=nbs; % number of boundary segments
    return
end
% E1 E2
d1=[
    0 0  % start parameter value
    1 1  % end parameter value
    1 3  % left hand region
    0 1 % right hand region
    ];
% E3 E4
d2=[
    0 0  % start parameter value
    1 1 % end parameter value
    2 1  % left hand region
    1 2 % right hand region
    ];
% E5 E6 E7 E8
d3=[
    0 0 0 0  % start parameter value
    1 1 1 1 % end parameter value
    2 1 0 0  % left hand region
    0 0 2 1 % right hand region
    ];
d4=[0 0 0
    1 1 1
    0 0 3
    3 3 0
    ];
d=[d1,d2,d3,d4];
bs1=bs(:)';

if find(bs1<1 | bs1>nbs)
    error(message('pde:lshapeg:InvalidBs'))
end

if nargin==1
    x=d(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2)
    error(message('pde:lshapeg:SizeBs'));
end
if ~isempty(s)
    
    % boundary segment 1
    ii=find(bs==1);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[X X],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[0 b],s(ii));
    end
    
    % boundary segment 2
    ii=find(bs==2);
    if ~isempty(ii)
        x(ii)=interp1([d(1,2),d(2,2)],[0 X],s(ii));
        y(ii)=interp1([d(1,2),d(2,2)],[b b],s(ii));
    end
    
    % boundary segment 3
    ii=find(bs==3);
    if ~isempty(ii)
        x(ii)=interp1([d(1,3),d(2,3)],[width width],s(ii));
        y(ii)=interp1([d(1,3),d(2,3)],[0 depthextend+faultdepth],s(ii));
    end
    
    % boundary segment 4
    ii=find(bs==4);
    if ~isempty(ii)
        x(ii)=interp1([d(1,4),d(2,4)],[0 width],s(ii));
        y(ii)=interp1([d(1,4),d(2,4)],[faultdepth+depthextend faultdepth+depthextend],s(ii));
    end
    % boundary segment 5
    ii=find(bs==5);
    if ~isempty(ii)
        x(ii)=interp1([d(1,5),d(2,5)],[0 width],s(ii));
        y(ii)=interp1([d(1,5),d(2,5)],[0 0],s(ii));
    end
    % boundary segment 6
    ii=find(bs==6);
    if ~isempty(ii)
        x(ii)=interp1([d(1,6),d(2,6)],[width X],s(ii));
        y(ii)=interp1([d(1,6),d(2,6)],[0 0],s(ii));
    end
    % boundary segment 7
    ii=find(bs==7);
    if ~isempty(ii)
        x(ii)=interp1([d(1,7),d(2,7)],[0 0],s(ii));
        y(ii)=interp1([d(1,7),d(2,7)],[0 depthextend+faultdepth],s(ii));
    end
    % boundary segment 8
    ii=find(bs==8);
    if ~isempty(ii)
        x(ii)=interp1([d(1,8),d(2,8)],[0 0],s(ii));
        y(ii)=interp1([d(1,8),d(2,8)],[depthextend+faultdepth b],s(ii));
    end
    % boundary segment 9
    ii=find(bs==9);
    if ~isempty(ii)
        x(ii)=interp1([d(1,9),d(2,9)],[0 0],s(ii));
        y(ii)=interp1([d(1,9),d(2,9)],[b Y],s(ii));
    end
    % boundary segment 10
    ii=find(bs==10);
    if ~isempty(ii)
        x(ii)=interp1([d(1,10),d(2,10)],[0 X],s(ii));
        y(ii)=interp1([d(1,10),d(2,10)],[Y Y],s(ii));
    end
    % boundary segment 11
    ii=find(bs==11);
    if ~isempty(ii)
        x(ii)=interp1([d(1,11),d(2,11)],[X X],s(ii));
        y(ii)=interp1([d(1,11),d(2,11)],[b Y],s(ii));
    end
end