function [x,y]=rectmutisubdomain(bs,s)
global mod faultdepth faultwidth
p = gcp;
noc = p.NumWorkers;
if noc==1
    error('This geometry can only run with Muti-workers')
elseif noc==2
    warning('More Workers please')
end
if faultdepth==0
    depthextend=100;
else
    depthextend=0;
end
if isempty(faultwidth)
    width=0.5;
else
    width=faultwidth;
end
X=mod.X;
Y=mod.Y;
nbs=2+3*noc;
if nargin==0
    x=nbs; % number of boundary segments
    return
end

d1=zeros(1,nbs); % start parameter value
d2=ones(1,nbs); % end parameter value
% d1(1,5)=1;
% d2(1,5)=0;
% d1(1,7:3:nbs-1)=1;
% d1(1,8:3:nbs)=1;
% d2(1,7:3:nbs-1)=0;
% d2(1,8:3:nbs)=0;

d3=zeros(1,nbs);
d4=zeros(1,nbs);
d3(1,1:5)=[2 1 1 1 1];% left hand region
d4(1,1:5)=[0 0 0 0 0];% right hand region
for i=1:noc-1
    d3(1,3+i*3:5+i*3)=(i+1);
    d4(1,3+i*3)=1;
    d4(1,4+i*3)=i+2;
    d4(1,5+i*3)=0;
end
d4(1,end-2:end-1)=1;
d=[d1;d2;d3;d4];
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

% point list
xl=zeros(1,3+2*noc);
yl=zeros(1,3+2*noc);
xl(1,1:5)=[0,width,X,X,0];
yl(1,1:5)=[0,0,0,Y,Y];
% divide fault zone in to noc-1 subdomains
dlength=faultdepth/(noc-1);
for i=1:noc-1
    xl(4+i*2)=width;
    xl(5+i*2)=0;
    yl(4+i*2:5+i*2)=dlength*i;
end
if ~isempty(s)
    % face 1
    % boundary segment 1
    for bsindex=1:4
        ii=find(bs==bsindex);
        if ~isempty(ii)
            x(ii)=interp1([d(1,1),d(2,1)],[xl(bsindex) xl(bsindex+1)],s(ii));
            y(ii)=interp1([d(1,1),d(2,1)],[yl(bsindex) yl(bsindex+1)],s(ii));
        end
    end

    ii=find(bs==5);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[xl(5) xl(end)],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[yl(5) yl(end)],s(ii));
    end
    % face 2
    ii=find(bs==6);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[xl(2) xl(6)],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[yl(2) yl(6)],s(ii));
    end
    ii=find(bs==7);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[xl(6) xl(7)],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[yl(6) yl(7)],s(ii));
    end
    ii=find(bs==8);
    if ~isempty(ii)
        x(ii)=interp1([d(1,1),d(2,1)],[xl(7) xl(1)],s(ii));
        y(ii)=interp1([d(1,1),d(2,1)],[yl(7) yl(1)],s(ii));
    end
    % face 3,4,5,.....
    for faceindex=3:noc
        ii=find(bs==faceindex*3);
        if ~isempty(ii)
            x(ii)=interp1([d(1,1),d(2,1)],[xl(faceindex*2) xl(faceindex*2+2)],s(ii));
            y(ii)=interp1([d(1,1),d(2,1)],[yl(faceindex*2) yl(faceindex*2+2)],s(ii));
        end
        ii=find(bs==faceindex*3+1);
        if ~isempty(ii)
            x(ii)=interp1([d(1,1),d(2,1)],[xl(faceindex*2+2) xl(faceindex*2+3)],s(ii));
            y(ii)=interp1([d(1,1),d(2,1)],[yl(faceindex*2+2) yl(faceindex*2+3)],s(ii));
        end
        ii=find(bs==faceindex*3+2);
        if ~isempty(ii)
            x(ii)=interp1([d(1,1),d(2,1)],[xl(faceindex*2+3) xl(faceindex*2+1)],s(ii));
            y(ii)=interp1([d(1,1),d(2,1)],[yl(faceindex*2+3) yl(faceindex*2+1)],s(ii));
        end
    end
end
end