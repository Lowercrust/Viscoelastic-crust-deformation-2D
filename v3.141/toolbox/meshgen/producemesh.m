% Producemesh using "generateMesh"
% BAD, not recommend
% mesh and boundary for finite element method
% input perameter 1 or 0; rectangular model geometry input 1, not rectangular input 0
function producemesh(model,varargin)
%% Generate Mesh
if nargin==2
    nn=varargin{1};
elseif nargin==1
    nn=500;% default setting
end
mesh=generateMesh(model,'Hmax',nn,'GeometricOrder','quadratic','Jiggle','on');
% if nargin==1
%     nn=varargin{1};
% elseif nargin==0
%     nn=7;% default setting
% elseif nargin>10
%     error('input number too large')
% else
%     error('inappropriate input argument')
% end
% for i=1:nn
%     [p,e,t] = refinemesh(g,p,e,t);
% end
end
