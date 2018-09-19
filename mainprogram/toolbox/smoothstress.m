% smooth stress 
function [stressyx,stressyz]=smoothstress(p,stressyx,stressyz,faultdepth0,faultdepth)

% smoothrange=30;
sr=50;
if faultdepth>faultdepth0 %&& faultdepth-faultdepth0>=1
    if faultdepth-faultdepth0<sr && faultdepth-faultdepth0>10
        smoothrange=(faultdepth-faultdepth0)*0.9;
    elseif faultdepth-faultdepth0<=10
        smoothrange=9;
    else
        smoothrange=sr;
    end
elseif faultdepth<faultdepth0
    if faultdepth0-faultdepth<sr*0.9
        if faultdepth0-faultdepth<10
            smoothrange=10;
        else
            smoothrange=(faultdepth0-faultdepth)*0.9;
        end
    elseif faultdepth0-faultdepth>sr*0.9
        smoothrange=sr*0.9;
    end
elseif faultdepth==faultdepth0
    faultdepth0=faultdepth-100;
    smoothrange=95;
end
nearfaultnode=find(p(1,:)<=smoothrange);
% disp(smoothrange);
nearfaulttipnode=find(abs(p(2,nearfaultnode)-(faultdepth0))<=smoothrange); % 10m
nearfaulttipnode=nearfaultnode(nearfaulttipnode);
x=p(1,nearfaulttipnode);
y=p(2,nearfaulttipnode);
xnodes = min(x):.1:max(x)+0.1;
ynodes = min(y):.1:max(y)+0.1;
zyx=stressyx(nearfaulttipnode);
zyz=stressyz(nearfaulttipnode);
zyx = RegularizeData3D(x,y,zyx,xnodes,ynodes,'smoothness',1);
zyz = RegularizeData3D(x,y,zyz,xnodes,ynodes,'smoothness',1);
[X,Y] = ndgrid(xnodes,ynodes);
interpzyx=griddedInterpolant(X,Y,zyx');
interpzyz=griddedInterpolant(X,Y,zyz');
zyx1=zeros(size(x,1),size(x,2));
zyz1=zeros(size(x,1),size(x,2));
for i=1:length(nearfaulttipnode)
    zyx1(i)=interpzyx(x(i),y(i));
    zyz1(i)=interpzyz(x(i),y(i));
end
stressyx(nearfaulttipnode)=zyx1';
stressyz(nearfaulttipnode)=zyz1';
end