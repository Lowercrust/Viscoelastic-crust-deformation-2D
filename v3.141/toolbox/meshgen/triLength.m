function tL=triLength(p,t)
sizet=length(t);
tL=zeros(sizet,4);
for i=1:sizet
    tL(i,1)=hypot(p(1,t(1,i))-p(1,t(2,i)),p(2,t(1,i))-p(2,t(2,i)));
    tL(i,2)=hypot(p(1,t(1,i))-p(1,t(3,i)),p(2,t(1,i))-p(2,t(3,i)));
    tL(i,3)=hypot(p(1,t(2,i))-p(1,t(3,i)),p(2,t(2,i))-p(2,t(3,i)));
    X=p(1,t(1:3,i));
    Y=p(2,t(1:3,i));
    tL(i,4)=polyarea(X,Y);
end