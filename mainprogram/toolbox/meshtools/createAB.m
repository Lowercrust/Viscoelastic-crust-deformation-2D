
function [A,B]=createAB(p,t)
np=size(p,2);
nt=size(t,2);
gpustatus(true)
A1=sparse(t(1:3,:),ones(3,1,'gpuArray')*(1:nt),1,np,nt);
A=gather(A1);
clear A1;
C=gpuArray(1./sum(A.'));
B1=sparse(1:np,1:np,C,np,np);
clear C
B=gather(B1);
clear B1
gpustatus(false)
end