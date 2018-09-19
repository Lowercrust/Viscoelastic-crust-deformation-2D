% create right hand side for mechanical equation
% input argument order 1: mesh, 2:dt, 3:stressyx, 4:stressyz, 5:esyxv,
% 6:esyzv;
function right=createright(msh,dt,varargin)
global modc fi
para = inputParser;
addRequired(para,'datatype');
addOptional(para,'tyx',[])
addOptional(para,'tyz',[])
addOptional(para,'esyxv',[])
addOptional(para,'esyzv',[])
addOptional(para,'esyxv2',[])
addOptional(para,'esyzv2',[])
parse(para,msh,dt,varargin{:})
[p,~,t] = meshToPet(msh);
if isempty(para.Results.esyxv)
    data=cell(2,1);
    data{2}=trigrad(p,t,para.Results.tyz,'output','dy','form','node');
    data{1}=trigrad(p,t,para.Results.tyx,'output','dx','form','node');
    data{4}=para.Results.esyzv2;
    data{3}=para.Results.esyxv2;
elseif isempty(para.Results.esyxv2)
    data=cell(4,1);
    data{2}=trigrad(p,t,para.Results.tyz,'output','dy','form','node');
    data{1}=trigrad(p,t,para.Results.tyx,'output','dx','form','node');
    data{4}=trigrad(p,t,para.Results.esyzv,'output','dy','form','node');
    data{3}=trigrad(p,t,para.Results.esyxv,'output','dx','form','node');
end
% stress
if length(modc.rhe)~=1
    p1i=fi{1};
    p2i=fi{2};
    right1=zeros(length(p),1);
    crhe=modc.rhe{1};
    mrhe=modc.rhe{2};
    dataplus=(data{1}+data{2});
    right1(p1i,1)=(dataplus(p1i,1))/(crhe.G*dt);
    right1(p2i,1)=(dataplus(p2i,1))/(mrhe.G*dt);
elseif length(modc.rhe)==1
    right1=(data{1}+data{2})/(modc.rhe.G*dt);
end
% shear strain
right2=data{3}+data{4};
right3=-right2+right1;
right=pdeintrpgpu(p,t,right3); % Interpolate from node data to triangle midpoint data
end

% end