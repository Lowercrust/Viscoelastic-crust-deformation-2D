% triangle surface area , optional input: point index, give the average surface area of triangle that connect to the point
function S=triareapt(p,t,varargin)
pc=cell(3,2);
if isempty(varargin)
    for i=1:3
        for j=1:2
            pc{j,i}=p(j,t(i,:));
        end
    end
    S=(pc{1,1}.*pc{2,2}+pc{1,2}.*pc{2,3}+pc{1,3}.*pc{2,1}-pc{1,3}.*pc{2,2}-pc{1,2}.*pc{2,1}-pc{1,1}.*pc{2,3}).*.5;
else
    [~,tindex]=find(t(1:3,:)==varargin{1});
    for i=1:3
        for j=1:2
            pc{j,i}=p(j,t(i,tindex));
        end
    end
    S=(pc{1,1}.*pc{2,2}+pc{1,2}.*pc{2,3}+pc{1,3}.*pc{2,1}-pc{1,3}.*pc{2,2}-pc{1,2}.*pc{2,1}-pc{1,1}.*pc{2,3}).*.5;
    S=mean(S);
end
end