function [Kc, Fc, null, ud, M] = applyBCNullSpacesa(K,A,F,Q,G,H,R,M)
% K = FEMatrices.K;
% A = FEMatrices.A;
% F = FEMatrices.F;
% Q = FEMatrices.Q;
% G = FEMatrices.G;
% H = FEMatrices.H;
% R = FEMatrices.R;
% M = FEMatrices.M;

[null,orth]=pdenullorth(H);
if size(orth,2)==0
    ud=zeros(size(K,2),1);
else
    ud=full(orth*((H*orth)\R));
end
Kc=K+A+Q;
Fc=null'*((F+G)-Kc*ud);
Kc=null'*Kc*null;
Kc = pde.PDEModel.checkForSymmetry(Kc);
M = null'*M*null;

end