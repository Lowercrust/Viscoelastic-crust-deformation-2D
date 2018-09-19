clear
close all
global mod
filename='anorthite(wet).xlsx';
[~,~,model]=xlsread('model.xlsx');
mod=cell2table(model(:,2)');
mod.Properties.VariableNames=model(:,1)';
mod=table2struct(mod);
mod.X=1;
mod.Y=1;
g=@rect;
model = createpde;
geometryFromEdges(model,g);
[p,e,t] = initmesh(g,'hmax',inf,'init','on'); 
for i=1:9
[p,e,t] = refinemesh(g,p,e,t);
end
model.Mesh= PettoMesh(p,e,t);
applyBoundaryCondition(model,'edge',3,'u',1);
applyBoundaryCondition(model,'edge',1,'u',0);
applyBoundaryCondition(model,'edge',2,'q',1);
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);
result = solvepde(model);
% tic
% FEM = assembleFEMatrices(model);
% toc

tic
[Kgpu,Fgpu]=createGlobalKFGPU(p,t,1,1);
toc
tic
[Kp,Fp]=createGlobalKFp(p,t,1,1);
toc
% time=cputime;
% [Ks,Fs]=createGlobalKFs(p,t,1,1);
% time=cputime-time;
% disp(time);
[Q,G,H,R]=assembleBoundary(model);
coefstruct = packCoefficientssa(model.EquationCoefficients);
A = formGlobalM2D(model, p, t, coefstruct,[],0,'a');
tic
[K, F] = formGlobalKF2D(model, p, t, coefstruct,[],0);
toc