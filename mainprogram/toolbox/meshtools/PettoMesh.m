function msh=PettoMesh(p,e,t)
tt = t';
tt(:,end) = [];
pt = p';
tr = triangulation(tt, p');
edges = tr.edges();
dx = pt(edges(:,2),1)-pt(edges(:,1),1);
dy = pt(edges(:,2),2)-pt(edges(:,1),2);
meshmaxsz = max(sqrt(sum(abs([dx dy]).^2,2)));
meshminsz = min(sqrt(sum(abs([dx dy]).^2,2)));
assoc = pde.FEMeshAssociation(t, e);
t(end,:) = [];
msh = pde.FEMesh(p, t, meshmaxsz, meshminsz, 'linear', assoc);
end