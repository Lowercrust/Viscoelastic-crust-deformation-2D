function L=edgeLength(p,e)
sizee=length(e);
L=zeros(sizee,1);
for i=1:sizee
    L(i)=hypot(p(1,e(1,i))-p(1,e(2,i)),p(2,e(1,i))-p(2,e(2,i)));
end
end