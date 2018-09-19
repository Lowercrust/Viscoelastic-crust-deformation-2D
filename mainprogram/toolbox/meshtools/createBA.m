function [A,B]=createBA(p,t)
np=size(p,2);
nt=size(t,2);
gpustatus(true)
A1=sparse(ones(3,1,'gpuArray')*(1:nt),t(1:3,:),1,nt,np);
A=gather(A1);
clear A1
gpustatus(false)
B=sparse(1:nt,1:nt,1./sum(A.'),nt,nt);
end