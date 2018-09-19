function [x,y]=rectwithtwosubdomain(bs,s)
global mod faultdepth
xs=10;
if faultdepth==0
    ys=100;
else
    ys=faultdepth;
end
X=mod.X;
Y=mod.Y;
nbs=10;
if nargin==0
    x=nbs; % number of boundary segments
    return
end
d=[%1 2 3 4 5 6 7 8 9 10
    0 0 0 0 0 0 0 0 0 0% start parameter value
    1 1 1 1 1 1 1 1 1 1 % end parameter value
    2 1 1 0 2 0 0 0 3 3% left hand region
    0 0 0 1 1 3 2 3 2 1% right hand region
    ];
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
        x(ii)=interp1([d(1,1),d(2,1)],[0 xs],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
    end
    % boundary segment 2
    ii=find(bs==2);
    if ~isempty(ii)
        x(ii)=interp1([d(1,2),d(2,2)],[xs X],s(ii));
        y(ii)=interp1([d(1,2),d(2,2)],[0 0],s(ii));
    end
    
    % boundary segment 3
    ii=find(bs==3);
    if ~isempty(ii)
        x(ii)=interp1([d(1,3),d(2,3)],[X X],s(ii));
        y(ii)=interp1([d(1,3),d(2,3)],[0 Y],s(ii));
    end
    
    % boundary segment 4
    ii=find(bs==4);
    if ~isempty(ii)
        x(ii)=interp1([d(1,4),d(2,4)],[xs X],s(ii));
        y(ii)=interp1([d(1,4),d(2,4)],[Y Y],s(ii));
    end
    
    % boundary segment 5
    ii=find(bs==5);
    if ~isempty(ii)
        x(ii)=interp1([d(1,5),d(2,5)],[xs xs],s(ii));
        y(ii)=interp1([d(1,5),d(2,5)],[0 ys],s(ii));
    end
    % boundary segment 6
    ii=find(bs==6);
    if ~isempty(ii)
        x(ii)=interp1([d(1,6),d(2,6)],[0 xs],s(ii));
        y(ii)=interp1([d(1,6),d(2,6)],[Y Y],s(ii));
    end
    % boundary segment 7
    ii=find(bs==7);
    if ~isempty(ii)
        x(ii)=interp1([d(1,7),d(2,7)],[0 0],s(ii));
        y(ii)=interp1([d(1,7),d(2,7)],[0 ys],s(ii));
    end
    % boundary segment 8
    ii=find(bs==8);
    if ~isempty(ii)
        x(ii)=interp1([d(1,8),d(2,8)],[0 0],s(ii));
        y(ii)=interp1([d(1,8),d(2,8)],[ys Y],s(ii));
    end
    % boundary segment 9
    ii=find(bs==9);
    if ~isempty(ii)
        x(ii)=interp1([d(1,9),d(2,9)],[0 xs],s(ii));
        y(ii)=interp1([d(1,9),d(2,9)],[ys ys],s(ii));
    end
    % boundary segment 10
    ii=find(bs==10);
    if ~isempty(ii)
        x(ii)=interp1([d(1,10),d(2,10)],[xs xs],s(ii));
        y(ii)=interp1([d(1,10),d(2,10)],[ys Y],s(ii));
    end
end