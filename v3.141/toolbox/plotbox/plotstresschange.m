function plotstresschange(stressbefore,stressafter,modeln,center)
% addpath()
global modc
figure('Position',[0 0 600 800]);
stresschange=stressbefore-stressafter;
% stresschange=smoothrange(modeln.Mesh.Nodes,double(stresschange),[0 modc.mod.CT-1000 modc.mod.X modc.mod.CT],50,0.001);
stresschange=smoothrange(modeln.Mesh.Nodes,double(stresschange),[0 0 modc.mod.X modc.mod.CT],50,0.00001);
stresschange=smoothrange(modeln.Mesh.Nodes,double(stresschange),[0 modc.mod.CT modc.mod.X modc.mod.Y],50,0.00001);
plus=stresschange;
minus=stresschange;
plus(plus<0)=0;
minus(minus>0)=0;
plus=log10(plus);
minus=-log10(abs(minus));
% center 
plus=plus-center;
minus=minus+center;

plus(plus<0)=0;
minus(minus>0)=0;
log=plus+minus;
% log(abs(log)<4)
pdeplot(modeln,'xydata',log) %,'contour','on','levels',10
colormap(brewermap([],'RdBu'))
% msc=max(abs(log));
msc=3;
caxis([-msc msc])
axis ij equal tight
title('+:stress drop -:stress increase')
% figure
% pdeplot(modeln,'xydata',stresschange) 
% colormap(brewermap([],'BrBG'))
% caxis([-10.^(msc-1) 10.^(msc-1)])
% axis ij
% xlim([0,100]);
% ylim([0,100]);
end

% smooth stress
function [out]=smoothrange(p,in,range,varargin)
node=findnodeinrect(p,range);
x=p(1,node);
y=p(2,node);
if isempty(varargin)
    xnodes = range(1):.1:range(3);
    ynodes = range(2):.1:range(4);
else
    xnodes = range(1):varargin{1}:range(3);
    ynodes = range(2):varargin{1}:range(4);
end
zyx1=in(node);
if length(varargin)~=2
    zyx = RegularizeData3D(x,y,zyx1,xnodes,ynodes,'smoothness',1);
else
    zyx = RegularizeData3D(x,y,zyx1,xnodes,ynodes,'smoothness',varargin{2},'interp','bicubic');
end
[X,Y] = ndgrid(xnodes,ynodes);
interpzyx=griddedInterpolant(X,Y,zyx');
zyx1=zeros(size(x,1),size(x,2));
for i=1:length(node)
    zyx1(i)=interpzyx(x(i),y(i));
end
out=in;
out(node)=zyx1';
end

