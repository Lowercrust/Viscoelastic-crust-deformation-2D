% xq2=0:25:mod.X;
% yq2=0:25:mod.Y;
% stressyxgrid=tri2grid(p,t,stressyx,xq2,yq2);
% tic;[Tfield,MaxdiffdT,MaxdiffdT0]=solvelinearTDPDE(Tfield,[0,dt],femodel);toc

tic;Fdisl=FPT(Pfield,Tfield,modelt.Mesh,'disl')';toc
% A = sprand(1000,1000,0.005);
% B = sprand(1000,1000,0.005);
% tic;Ag= gpuArray(A);toc
% Bg= gpuArray(B);
% timeit(@() A\B)
% timeit(@() Ag\Bg)