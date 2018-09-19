function un=pdeprtnigpu(p,t,ut)
global AB
%PDEPRTNI Interpolate from triangle data midpoints to node data.
%
%       UN=PDEPRTNI(P,T,UT) gives linearly interpolated values
%       at node point from the values at triangle mid points.
%
%       The geometry of the PDE problem is given by the triangle data P
%       and T. Details under INITMESH.
%
%       Let N be the dimension of the PDE system, and NP the number of
%       node points and NT the number of triangles. The components
%       of triangle data in UT are stored as N rows of length NT.
%       The components of the node data are stored in UN as N columns
%       of length NP.
%
%       See also ASSEMPDE, INITMESH, PDEINTRP

%       A. Nordmark 12-14-94.
%       Copyright 1994-2003 The MathWorks, Inc.
if isempty(AB)
    np=size(p,2);
    nt=size(t,2);
    if size(ut,2) ~= nt
        error(message('pde:pdeprtni:NumColsUt'));
    end
    % waitGPU;
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
    un=B*A*(ut.');
else
    un=AB*(ut.');
end
