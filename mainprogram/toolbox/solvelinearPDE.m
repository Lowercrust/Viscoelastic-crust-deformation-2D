% input argument order
% right(double matrix), boundary (struct),  K(double sparse)
function v=solvelinearPDE(model,varargin)
[p,~,t]=meshToPet(model.Mesh);
c=model.EquationCoefficients.CoefficientAssignments.c;
% coefstruct = packCoefficientssa(model.EquationCoefficients); % 2.5s
A=sparse(size(p,2),size(p,2),0);
M=sparse(size(p,2),size(p,2),0);
% A = formGlobalM2D(model, p, t, coefstruct,[],0,'a'); % 3.5s
% M = formGlobalM2D(model, p, t, coefstruct,[],0,'m'); % 3.5s
switch nargin
    case 1
        f=model.EquationCoefficients.CoefficientAssignments.f;
        [K,F]=createGlobalKFGPU(p,t,'K',c,'F',f); % 6s;
        bmatrix=struct('Q',{},'G',{},'H',{},'R',{});
        [bmatrix(1).Q,bmatrix(1).G,bmatrix(1).H,bmatrix(1).R]=assembleBoundary(model);
    case 2 % 
        fmatrix=varargin{1};
        [K,F]=createGlobalKFGPU(p,t,'K',c,'F',fmatrix);
        bmatrix=struct('Q',{},'G',{},'H',{},'R',{});
        [bmatrix(1).Q,bmatrix(1).G,bmatrix(1).H,bmatrix(1).R]=assembleBoundary(model);
    case 3
        fmatrix=varargin{1};
        bmatrix=varargin{2};
        [K,F]=createGlobalKFGPU(p,t,'K',c,'F',fmatrix);
    case 4 % not cal K
        fmatrix=varargin{1};
        bmatrix=varargin{2};
        K=varargin{3};
        F=createGlobalKFGPU(p,t,'F',fmatrix);
end

[Kc, Fc, B, ud, ~] = applyBCNullSpacesa(K,A,F,bmatrix.Q,bmatrix.G,bmatrix.H,bmatrix.R,M); % 0.7s
v= B*(Kc\Fc) +ud; % 1.0s
% varlist = {'Kc','Fc'};
% clear(varlist{:})
% v = gather(v);
% varlist = {'Kc','Fc','B','ud'};
% clear(varlist{:})
% result= createPDEResults(model,v);
end



