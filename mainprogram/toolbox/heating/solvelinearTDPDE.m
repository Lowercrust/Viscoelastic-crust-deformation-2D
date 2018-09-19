function [Tfield,varargout]=solvelinearTDPDE(u0,tlist,varargin)
global femodel Tfield0
femodel=varargin{1}; % femodel is global in this function
% [p,e,t] = meshToPet(model.Mesh);
% u0 = imposeBConICfirstOrderODE(self,u0,tlist);
nu=size(u0,2);
%Following pragmas, %#function, are needed for deployment.
%#function pde.internal.odefunchandles.firstOrderODEf
fcn  = str2func('ODEf');
%#function pde.internal.odefunchandles.firstOrderODEdf
dfcn = str2func('ODEdf');
%#function pde.internal.odefunchandles.firstOrderODEm
mfcn = str2func('ODEm');
tsecondOrder=false;
odeoptions = constructODEoptionssa(u0,1e-3,1e-6,'off',nu,dfcn,mfcn,tsecondOrder);
uu0=imposeBConICfirstOrderODEsa(u0)';
[~,uu]=ode15s_test(fcn,tlist,uu0,odeoptions);
numCompTimes = size(uu, 1);
DoFIndex=1:size(uu,2); % vector of 1 to number of degree of freedom

if(length(tlist)==2)
    % ODE15S behaves differently if length(tlist)==2
    numCompTimes = 2;
    u1(DoFIndex,:)    = uu([1 size(uu,1)],DoFIndex)';
else
    u1(DoFIndex,:)    = uu(:,DoFIndex)';
end

% u1(DoFIndex,:)    = uu([1 size(uu,1)],DoFIndex)';
u1=femodel.B*u1+femodel.ud*ones(1,numCompTimes);
Tfield=u1(:,end)';
dT=u1(:,end)'-u1(:,1)';
dT0=u1(:,end)'-Tfield0;
[MaxTdiff,MI]=max(abs(dT));
MaxTdiff=MaxTdiff*sign(dT(MI));
[MaxTdiff0,MI]=max(abs(dT0));
MaxTdiff0=MaxTdiff0*sign(dT0(MI));
varargout{1}=MaxTdiff;
varargout{2}=MaxTdiff0;
end

function f=ODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015 The MathWorks, Inc.
global femodel
gpustatus(true)
Kg=gpuArray(femodel.K);
u=gpuArray(u);
F=gpuArray(femodel.F);
f1=-Kg*u+F;
clear Kg F u
f=gather(f1);
clear f1
gpustatus(false)
end

function df=ODEdf(t,u)
%firstOrderODEdf - Jacobian function for first-order ODE system
%Computes jacobian of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015 The MathWorks, Inc.

global femodel
df=-femodel.K;
end



function m=ODEm(t,u)
%firstOrderODEm - Mass matrix function for first-order ODE system
%Computes the mass matrix of discretized PDE with first-order time
%derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015 The MathWorks, Inc.
global femodel

% m=femodel.B'*femodel.Mass*femodel.B;
m=femodel.m;
end


function odeoptions = constructODEoptionssa(uu0,rtol,atol,stats,nu,dfcn,mfcn,tsecondOrder)
odeoptions = odeset;
% If doJacFiniteDiff is true, ode15s computes the Jacobian by finite
% difference. This is useful for testing purposes.
doJacFiniteDiff=false;
if(doJacFiniteDiff)
    jac0=fcn(0,uu0);
    [ir,ic]=find(jac0);
    jpat=sparse(ir,ic,ones(size(ir,1),1), size(jac0,1), size(jac0,2));
    odeoptions=odeset(odeoptions, 'JPattern', jpat);
else
    odeoptions=odeset(odeoptions, 'Jacobian',dfcn);
end
odeoptions=odeset(odeoptions,'Mass',mfcn);
odeoptions=odeset(odeoptions,'RelTol',rtol);
if (tsecondOrder)
    % Mask out dudt
    atol=atol*ones(2*nu,1);
    atol(nu+1:nu+nu)=Inf*ones(nu,1);
    odeoptions=odeset(odeoptions,'MaxOrder',2);
end

odeoptions=odeset(odeoptions,'AbsTol',atol);
odeoptions=odeset(odeoptions,'Stats',stats);
odeoptions=odeset(odeoptions,'MassSingular','no');
end

function u0 = imposeBConICfirstOrderODEsa(u0,B)
global femodel
u0=femodel.B'*u0';
end