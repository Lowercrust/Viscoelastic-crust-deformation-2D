% only for rectangle geometry and subdomain
% create t with domainindex
function t=domainindex(p,e)
delTri = pdeDelTri(p');
t = delTri.getTriangles();
S=triareapt(p,t);
t0=t;
t(:,S<=0.001)=[];% remove triangle with small surface area
minS=min(S);
if (size(t0,2)-size(t,2))~=0
    disp([num2str(size(t0,2)-size(t,2)) ' small triangle (<0.001) has been removed'])
    disp(['minum triangle surface area before:' num2str(minS)])
    S(:,S<=0.001)=[];
    disp(['minum triangle surface area after:' num2str(min(S))])
end
% if isempty(varargin)
% center of each tri.
rt=reshape(t(1:3,:),[1,size(t,2)*3]);
tx=p(1,rt);
ty=p(2,rt);
tx=reshape(tx,[3,size(t,2)]);
ty=reshape(ty,[3,size(t,2)]);
cx=mean(tx,1);
cy=mean(ty,1);
t(4,:)=1;
for i=2:max(max(e(6:7,:)))
    index=find(e(6,:)==i);
    ex=[p(1,e(1,index)),p(1,e(2,index))];
    ey=[p(2,e(1,index)),p(2,e(2,index))];
    maxx=max(ex);
    maxy=max(ey);
    minx=min(ex);
    miny=min(ey);
    %         disp([size(cx),i,ex])
    if ~isempty(ex)
        indextx= find(cx<maxx & cx>minx);
        indexty= find(cy(indextx)<maxy & cy(indextx)>miny);
        t(4,indextx(indexty))=i;
    end
end
% else
%     t(4,:)=varargin{1};
% end
end

