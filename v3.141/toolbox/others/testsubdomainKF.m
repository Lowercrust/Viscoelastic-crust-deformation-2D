% 
index1=find(l1.t(4,:)==1);
point1=l1.t(1:3,index1);
point1=reshape(point1,length(point1)*3,1);
point1=unique(point1);
pface1=l1.p(:,point1);
eface1=createedgematrix(g,pface1);
tface1=domainindex(pface1,eface1);
index1d=find(tface1(4,:)~=1);
tface1(:,index1d)=[];
coefstructface1=l1.coefstruct;
coefstructface1.ElementsInSubdomain=(1:length(tface1))';
coefstructface1.NumElementsInSubdomain=length(tface1);
% face1mesh=PettoMesh(pface1,eface1,tface1);
model.Mesh=PettoMesh(pface1,eface1,tface1);
[K1, F1] = formGlobalKF2D(model, pface1, tface1, coefstructface1,[],0);
model.Mesh=msh;
[K, F] = formGlobalKF2D(l1.thePde, l1.p, l1.t, l1.coefstruct,[],0);

