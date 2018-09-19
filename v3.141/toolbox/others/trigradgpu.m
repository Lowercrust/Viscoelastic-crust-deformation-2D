% evaluate gradient for a solution
function varargout=trigradgpu(p,t,u,varargin)
tic
x1=p(1,t(1,:));
x2=p(1,t(2,:));
x3=p(1,t(3,:));
y1=p(2,t(1,:));
y2=p(2,t(2,:));
y3=p(2,t(3,:));
u1=u(t(1,:));
u2=u(t(2,:));
u3=u(t(3,:));
A=gpuArray([x2-x1;y2-y1;(u2-u1)']);
B=gpuArray([x3-x1;y3-y1;(u3-u1)']); 
C=cross(A,B,1); % Normal vector of each triangle
C=gather(C);
toc
if nargout == 2
    gradx=C(1,:)./C(3,:); % gradient of each triangle
    grady=C(2,:)./C(3,:);
    varargout{1}=-pdeprtni(p,t,gradx); % interp. tri-data to node data.
    varargout{2}=-pdeprtni(p,t,grady);
elseif nargout == 1
    if strcmp(varargin{1},'dx')
        gradx=C(1,:)./C(3,:); % gradient of each triangle
        varargout{1}=-pdeprtni(p,t,gradx); % interp. tri-data to node data.
    elseif strcmp(varargin{1},'dy')
        grady=C(2,:)./C(3,:);
        varargout{1}=-pdeprtni(p,t,grady);
    end
end
varlist = {'A','B','C','gradx','grady'};
clear(varlist{:})
end