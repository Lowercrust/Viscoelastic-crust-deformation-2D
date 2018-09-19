% refine a mesh around a point [x0,y0] using matlab function "refinemesh"
function msh=pointrefine(g,msh,x0,y0,r,target) %p1,e1,t1
disp(msh.MinElementSize)
refinecount=0;
while msh.MinElementSize>target
    [p,e,t]=meshToPet(msh);
    sizep=length(p);
    % msh=PettoMesh(p,e,t);
    index=zeros(sizep,1);
    j=0;
    % find meshs to refine
    for i=1:sizep
        if hypot((p(1,i)-x0),(p(2,i)-y0))<r % p(1,i)<=x1 && p(1,i)>=x0 && p(2,i)<=y1 && p(2,i)>=y0
            j=j+1;
            index(j,1)=i;
        end
    end
    if r>target
        if r>msh.MaxElementSize;
            r=msh.MaxElementSize;
        else
            r=r/1.5;
        end
    end
    index(j+1:end)=[];
    
    sizet=size(t,2);
    j=0;
    it=zeros(sizet,1);
    %disp(t)
    for i=1:sizet
        tri=t(1:3,i);
        if ~isempty(find(index==tri(1),1)) && ~isempty(find(index==tri(2),1)) && ~isempty(find(index==tri(3),1))
            j=j+1;
            it(j,1)=i;
        end
        
    end
%     disp(j)
    it(j+1:end)=[];
    tic
    [p,e,t] = refinemesh(g,p,e,t,it);
    refinecount=refinecount+1;
    MinElementSize0=msh.MinElementSize;
    msh=PettoMesh(p,e,t);
    et=toc;
    meshsize=size(msh.Nodes,2);
    disp(['time ' num2str(et) ' Num of Nodes' num2str(meshsize) ' MinElementsize ' num2str(msh.MinElementSize)])

    if MinElementSize0==msh.MinElementSize || refinecount>10
        %         pdemesh(p,e,t)
        %         hold on
        %         scatter([x0n,x1n],[y0n,y1n])
        %         hold off
        pause
        break
    end
end
end