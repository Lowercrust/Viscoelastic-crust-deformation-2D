function varargout=createGlobalKFGPU(p,t,varargin)
pc=cell(3,2);
for i=1:3
    for j=1:2
        pc{j,i}=p(j,t(i,:));
    end
end
% center of triangles
xpts = (pc{1,1} + pc{1,2} + pc{1,3})/3;
ypts = (pc{2,1} + pc{2,2} + pc{2,3})/3;
% Area of triangles
S=(pc{1,1}.*pc{2,2}+pc{1,2}.*pc{2,3}+pc{1,3}.*pc{2,1}-pc{1,3}.*pc{2,2}-pc{1,2}.*pc{2,1}-pc{1,1}.*pc{2,3}).*.5;
% center of each triangle
region.x=xpts;
region.y=ypts;
varargout=cell(2);
for l=1:2:length(varargin)
    
    c=varargin{l+1};
    if strcmp(varargin{l},'K')
        if isa(c,'double')
            c1=ones(1,size(t,2)).*c;
        elseif isa(c,'function_handle')
            c1=c(region);
        elseif isa(c,'scatteredInterpolant')
            c1=c(xpts,ypts);
        end
        %% K matrix
        % index of Global K ;
        % edge vector of triangles
        vx=zeros(3,size(t,2));
        vy=zeros(3,size(t,2));
        vcell=cell(1,9);
        vx(1,:)=pc{1,3}-pc{1,2};
        vx(2,:)=pc{1,1}-pc{1,3};
        vx(3,:)=pc{1,2}-pc{1,1};
        vy(1,:)=pc{2,3}-pc{2,2};
        vy(2,:)=pc{2,1}-pc{2,3};
        vy(3,:)=pc{2,2}-pc{2,1};
        indexcellx=cell(1,9);
        indexcelly=cell(1,9);
        for i=1:3
            for j=1:3
                indexcellx{(i-1)*3+j}=t(i,:);
                indexcelly{(i-1)*3+j}=t(j,:);
                vcell{(i-1)*3+j}=(vx(i,:).*vx(j,:)+ vy(i,:).*vy(j,:)).*c1./(4*S);
            end
        end
        indexx=uint32(cell2mat(indexcellx));
        indexy=uint32(cell2mat(indexcelly));
        v=cell2mat(vcell);
        gpustatus(true)
        indexx=gpuArray(indexx');
        indexy=gpuArray(indexy');
        v=gpuArray(v');
        K=sparse(indexx,indexy,v);
        vars={'indexx','indexy','v'};
        clear(vars{:});
        K1=gather(K);
        clear K
        gpustatus(false)
        varargout{l}=K1;
    elseif strcmp(varargin{l},'F')
        f=varargin{l+1};
        if isa(f,'double')
            if length(f)==1
                f1=ones(1,size(t,2)).*f;
            elseif length(f)==size(t,2)
                f1=f;
            end
        elseif isa(f,'function_handle')
            f1=f(region);
        elseif isa(f,'scatteredInterpolant')
            f1=f(xpts,ypts);
        end
        %% F matrix
        fcell=cell(1,3);
        % index of Global F;
        indexcellf=cell(1,3);
        for i=1:3
            indexcellf{i}=t(i,:);
            fcell{i}=f1.*S/3;
        end
        indexf=single(cell2mat(indexcellf));
        ff=cell2mat(fcell);
        %         datasize=totalBytes(indexf)+totalBytes(ff);
        %         waitGPU;
        gpustatus(true)
        indexf=gpuArray(indexf);
        ff=gpuArray(ff);
        F=sparse(indexf',gpuArray(ones(length(indexf),1)),ff);
        vars={'indexf','ff'};
        clear(vars{:});
        F1=gather(F);
        clear F
        gpustatus(false)
        F=full(F1);
        varargout{(l+1)/2}=F;
    end
end
%%
varlist = {'indexx','indexy','v','indexf','ff'};
clear(varlist{:})
end