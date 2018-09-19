function [x,y]=crustmantletwoface(bs,s)
global modc faultdepth
X=modc.mod.X;
Y=modc.mod.Y;
b=modc.mod.CT;
if isempty(faultdepth)
    faultdepth=30000;
end
nbs=8;
if nargin==0,
    x=nbs; % number of boundary segments
    return
end

d=[
    0 0 0 0 0 0 0 0 % start parameter value
    1 1 1 1 1 1 1 1% end parameter value
    0 0 2 1 1 0 0 2% left hand region
    1 1 1 0 0 2 2 0% right hand region
    ];

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
if m==1 && n==1,
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
    error(message('pde:lshapeg:SizeBs'));
end
if ~isempty(s),
    
    % boundary segment 1
    ii=find(bs==1);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[0 faultdepth],s(ii));
    end
    
    % boundary segment 2
    ii=find(bs==2);
    if ~isempty(ii)
        x(ii)=interp1([d(1,2),d(2,2)],[0 0],s(ii));
        y(ii)=interp1([d(1,2),d(2,2)],[faultdepth b],s(ii));
    end
    
    % boundary segment 3
    ii=find(bs==3);
    if ~isempty(ii)
        x(ii)=interp1([d(1,3),d(2,3)],[0 X],s(ii));
        y(ii)=interp1([d(1,3),d(2,3)],[b b],s(ii));
    end
    
    % boundary segment 4
    ii=find(bs==4);
    if ~isempty(ii)
        x(ii)=interp1([d(1,4),d(2,4)],[X X],s(ii));
        y(ii)=interp1([d(1,4),d(2,4)],[0 b],s(ii));
    end
    
    % boundary segment 5
    ii=find(bs==5);
    if ~isempty(ii)
        x(ii)=interp1([d(1,5),d(2,5)],[0 X],s(ii));
        y(ii)=interp1([d(1,5),d(2,5)],[0 0],s(ii));
    end
    % boundary segment 6
    ii=find(bs==6);
    if ~isempty(ii)
        x(ii)=interp1([d(1,6),d(2,6)],[0 0],s(ii));
        y(ii)=interp1([d(1,6),d(2,6)],[b Y],s(ii));
    end
    % boundary segment 7
    ii=find(bs==7);
    if ~isempty(ii)
        x(ii)=interp1([d(1,7),d(2,7)],[0 X],s(ii));
        y(ii)=interp1([d(1,7),d(2,7)],[Y Y],s(ii));
    end
    % boundary segment 8
    ii=find(bs==8);
    if ~isempty(ii)
        x(ii)=interp1([d(1,8),d(2,8)],[X X],s(ii));
        y(ii)=interp1([d(1,8),d(2,8)],[b Y],s(ii));
    end
end