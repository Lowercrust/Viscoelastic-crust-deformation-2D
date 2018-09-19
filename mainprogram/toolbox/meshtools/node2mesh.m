% convert nodes data to mesh data
% input p: nodes and g, geometry
function msh=node2mesh(p,g)
% create edge
e=createedgematrix(g,p);
% create triangle
t=domainindex(p,e);
msh=PettoMesh(p,e,t);
end