function [pindex,x1,y1]=grid2pindex(x,y,p)
pindex=zeros(length(x),length(y));
x1=pindex;
y1=pindex;
for i=1:length(x)
    for j=1:length(y)
        [~,mini]=min((p(1,:)-x(i)).^2+(p(2,:)-y(j)).^2);
        pindex(i,j)=mini;
        x1(i,j)=p(1,mini);
        y1(i,j)=p(2,mini);
    end
    disp(i)
end
end