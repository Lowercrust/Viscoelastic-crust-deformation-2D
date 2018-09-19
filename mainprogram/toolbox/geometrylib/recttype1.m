function [x,y]=recttype1(bs,s)
global mod pointlist faultwidth
xs=faultwidth;
ys=pointlist(end-1);
ys0=pointlist(2);
X=mod.X;
Y=mod.Y;
nbs=9;
if nargin==0
    x=nbs; % number of boundary segments
    return
end
% E1 E2 E3
d1=[
    0 0 0 % start parameter value
    1 1 1 % end parameter value
    1 1 0 % left ys0and region
    0 0 1 % rigys0t ys0and region
    ];
% E4 E5 E6
d2=[
    0 0 0 % start parameter value
    1 1 1 % end parameter value
    2 2 1 % left ys0and region
    1 1 2 % rigys0t ys0and region
    ];
% E7 E8 E9
d3=[
    0 0 0  % start parameter value
    1 1 1  % end parameter value
    0 2 1  % left ys0and region
    1 0 0  % rigys0t ys0and region
    ];

d=[d1,d2,d3];
bs1=bs(:)';

if find(bs1<1 | bs1>nbs)
    error(message('pde:lsys0apeg:InvalidBs'))
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
    error(message('pde:lsys0apeg:SizeBs'));
end
if ~isempty(s)
    
    % boundary segment 1
    ii=find(bs==1);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[0 X],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[0 0],s(ii));
    end
    
    % boundary segment 2
    ii=find(bs==2);
    if ~isempty(ii)
        x(ii)=interp1([d(1,2),d(2,2)],[X X],s(ii));
        y(ii)=interp1([d(1,2),d(2,2)],[0 Y],s(ii));
    end
    
    % boundary segment 3
    ii=find(bs==3);
    if ~isempty(ii)
        x(ii)=interp1([d(1,3),d(2,3)],[0 X],s(ii));
        y(ii)=interp1([d(1,3),d(2,3)],[Y Y],s(ii));
    end
    
    % boundary segment 4
    ii=find(bs==4);
    if ~isempty(ii)
        x(ii)=interp1([d(1,4),d(2,4)],[0 xs],s(ii));
        y(ii)=interp1([d(1,4),d(2,4)],[ys ys],s(ii));
    end
    % boundary segment 5
    ii=find(bs==5);
    if ~isempty(ii)
        x(ii)=interp1([d(1,5),d(2,5)],[xs xs],s(ii));
        y(ii)=interp1([d(1,5),d(2,5)],[ys ys+ys0],s(ii));
    end
    % boundary segment 6
    ii=find(bs==6);
    if ~isempty(ii)
        x(ii)=interp1([d(1,6),d(2,6)],[0 xs],s(ii));
        y(ii)=interp1([d(1,6),d(2,6)],[ys+ys0 ys+ys0],s(ii));
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
        y(ii)=interp1([d(1,8),d(2,8)],[ys ys+ys0],s(ii));
    end
    % boundary segment 9
    ii=find(bs==9);
    if ~isempty(ii)
        x(ii)=interp1([d(1,9),d(2,9)],[0 0],s(ii));
        y(ii)=interp1([d(1,9),d(2,9)],[ys+ys0 Y],s(ii));
    end
    
end
