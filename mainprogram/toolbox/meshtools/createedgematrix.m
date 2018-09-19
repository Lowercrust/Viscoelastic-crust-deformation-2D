% create edge matrix using model geometry and node coordinate data
function edge=createedgematrix(g,p1)
ne=g();
[x0,y0]=g(1:ne,zeros(1,ne)); % start points;
[x1,y1]=g(1:ne,ones(1,ne)); % end points;

edge=cell(1,ne);
domaininfo=g(1:ne);
p1single=single(p1);
x0=single(x0);
y0=single(y0);
x1=single(x1);
y1=single(y1);
for i=1:ne
    % find points indice of edge i
    if x0(i)==x1(i)
        [~,col0]=find(p1single(1,:)==x0(i));
        [~,col1]=find(p1single(2,col0)<=max(y0(i),y1(i))+1e-3 & p1single(2,col0)>=min(y0(i),y1(i))-1e-3);
        [arclength,col0]=sortdata(p1single(2,col0(col1)),col0(col1)); % sort col0
        arclength=arclength-min(y0(i),y1(i));
        arclength=arclength/abs(y1(i)-y0(i));
    else
        [~,col0]=find(p1single(2,:)==y0(i));
        [~,col1]=find(p1single(1,col0)<=max(x0(i),x1(i))+1e-3 & p1single(1,col0)>=min(x0(i),x1(i))-1e-3);
        [arclength,col0]=sortdata(p1single(1,col0(col1)),col0(col1)); % sort col0
        arclength=arclength-min(x0(i),x1(i));
        arclength=arclength/abs(x1(i)-x0(i));
    end
    edge0=zeros(7,size(col1,2)-1);
    if size(col1,2)==2
        edge0(3,:)=0;
        edge0(4,:)=1;
    elseif size(col1,2)>2
        edge0(3,:)=arclength(1:end-1);
        edge0(4,:)=arclength(2:end);
    end
    edge0(1,:)=col0(1:end-1);
    edge0(2,:)=col0(2:end);
    edge0(5,:)=ones(1,size(col1,2)-1)*i;
    edge0(6,:)=domaininfo(3,i)*ones(1,size(col1,2)-1);
    edge0(7,:)=domaininfo(4,i)*ones(1,size(col1,2)-1);
    edge{i}=edge0;    
end
edge=cell2mat(edge);