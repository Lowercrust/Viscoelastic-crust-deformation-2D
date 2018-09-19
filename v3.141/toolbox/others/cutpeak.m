function stress=cutpeak(model,stress,x0,x1,y0,y1)
% global faultdepth
index=zeros(length(stress),1);
j=0;
for i=1:length(stress)
    if model.Mesh.Nodes(2,i)<y1 && model.Mesh.Nodes(2,i)>y0 && model.Mesh.Nodes(1,i)<x1 && model.Mesh.Nodes(1,i)>x0
        j=j+1;
        index(j)=i;
    end
end
% disp(j)
index(j+1:end,:)=[];
% index=find(model.Mesh.Nodes(2,:)<faultdepth);
[~,~,t]=meshToPet(model.Mesh);
cut=zeros(length(stress),1);
for i=1:length(index)
    [row,col]=find(t(1:3,:)==index(i));
    adpindex=zeros(2,length(col));
    % find adjecent points;
    for j=1:length(col)
        tt=t(1:3,col(j));
        tt(row(j))=[];
        adpindex(:,j)=tt;
    end
    adpindex = reshape(adpindex,[length(col)*2,1]);
    meanstress=mean(stress(adpindex));
    evaluatestress=stress(index(i));
    if abs(evaluatestress-meanstress)>1e4
                cut(index(i))=evaluatestress-meanstress;
%         stress(index(i))=meanstress;
    end
end
stress=stress-cut;
end