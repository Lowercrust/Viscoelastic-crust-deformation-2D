% evaluate gradient on the center of each triangle
% output: node data
function varargout=trigrad(p,t,u,varargin)
% ni=nargin;
para = inputParser;
defaultoutput='dxdy';
defaultform='tri';
addRequired(para,'p',@isnumeric);
addRequired(para,'t',@isnumeric);
addRequired(para,'u',@isnumeric);
addOptional(para,'output',defaultoutput)
addOptional(para,'form',defaultform,@ischar)
parse(para,p,t,u,varargin{:})
% waitGPU;
gpustatus(true);
% p=gpuArray(p);
% u=gpuArray(u);
% A=[p(1,t(2,:))-p(1,t(1,:));p(2,t(2,:))-p(2,t(1,:));(u(t(2,:))-u(t(1,:)))'];
x1=gpuArray(p(1,t(1,:)));
x2=gpuArray(p(1,t(2,:)));
x3=gpuArray(p(1,t(3,:)));
y1=gpuArray(p(2,t(1,:)));
y2=gpuArray(p(2,t(2,:)));
y3=gpuArray(p(2,t(3,:)));
u1=gpuArray(u(t(1,:)));
u2=gpuArray(u(t(2,:)));
u3=gpuArray(u(t(3,:)));
A=[x2-x1;y2-y1;(u2-u1)'];
B=[x3-x1;y3-y1;(u3-u1)'];
varlist={'x1','x2','x3','y1','y2','y3','u1','u2','u3',};
clear(varlist{:});
C1=cross(A,B,1);
C=gather(C1); % Normal vector of each triangle
varlist={'A','B','C1'};
clear(varlist{:});
gpustatus(false)
switch para.Results.output
    case 'dx'
        gradx=C(1,:)./C(3,:); % gradient of each triangle
        varargout{1}=gradx;
    case 'dy'
        grady=C(2,:)./C(3,:);
        varargout{1}=grady;
    case 'dxdy'
        gradx=C(1,:)./C(3,:); % gradient of each triangle
        grady=C(2,:)./C(3,:);
        varargout{1}=gradx;
        varargout{2}=grady;
end
if strcmp(para.Results.form,'node')
    for i=1:length(varargout)
        varargout{i}=-pdeprtnigpu(p,t,varargout{i});
    end
end
% disp(nargout)
% if ni==0 % default setting: output x and y gradiant on tri-center;
%     gradx=C(1,:)./C(3,:); % gradient of each triangle
%     grady=C(2,:)./C(3,:);
%     varargout{1}=gradx;
%     varargout{2}=grady;
% elseif ni~=0

% if strcmp(varargin{1},'dx')
%     gradx=C(1,:)./C(3,:); % gradient of each triangle
%     varargout{1}=gradx;
%     varargout{1}=-pdeprtni(p,t,gradx); % interp. tri-data to node data.
% elseif strcmp(varargin{1},'dy')
%     grady=C(2,:)./C(3,:);
%     varargout{1}=grady;
%     varargout{1}=-pdeprtni(p,t,grady);
% end
% end

end