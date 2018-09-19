
function  [stressyx,stressyz,u]=stressdrop(stressyx,stressyz,modeln,gn,varargin)
global modc faultdepth fi
%time=cputime;
% calculate the solutions on a finer mesh
%% find fault depth on corser mesh
% msh=model.Mesh;
%reflesh modeln
% modeln.Mesh=msh;
if ~isempty(varargin{1})
    eqtype=varargin{1};
    faultdepth=findfaultdepth(modeln.Mesh,stressyx,eqtype);
else
    faultdepth=findfaultdepth(modeln.Mesh,stressyx);
end

modeln.Mesh=node2mesh(modc.Nodes,gn);
% p1=modeln.Mesh.Nodes;
% e1=createedgematrix(gn,p1);
% t1=domainindex(p1,e1);
% modeln.Mesh=PettoMesh(p1,e1,t1);
% disp(faultdepth);
%% refinemesh around the faultdepth
% eqmsh=generateMesh2d(modeln,gn,'pointrefinement');
% modeln.Mesh=eqmsh;
u=solvelinearPDE(modeln);
[p,~,t] = meshToPet(modeln.Mesh);
[esyxe,esyze]=trigrad(p,t,u,'form','node'); % elastic shear strain
stressdropyx=modc.rhe{1}.G*esyxe; % stressdrop after the earthquake yx
stressdropyz=modc.rhe{1}.G*esyze; % yz
if size(modc.rhe,1)==2
    stressdropyx(fi{2})=modc.rhe{2}.G*esyxe(fi{2}); % stressdrop after the earthquake yx
    stressdropyz(fi{2})=modc.rhe{2}.G*esyze(fi{2}); % yz
end
stressyx=stressyx-stressdropyx;
stressyz=stressyz-stressdropyz;
end

