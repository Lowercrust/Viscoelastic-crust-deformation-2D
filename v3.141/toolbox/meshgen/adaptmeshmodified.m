function [u,p,e,t,code]=adaptmeshmodified(g,b,c,a,f,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24)
%ADAPTMESH Adaptive mesh generation and PDE solution.
%
%       [U,P,E,T]=ADAPTMESH(G,B,C,A,F,P1,V1,...) performs adaptive mesh
%       generation and PDE solution.  The large number of possible input
%       option is handled using property-value pair arguments.  The first
%       five arguments G, B, C, A, and F are not optional.
%
%       The function produces a solution U to the elliptic scalar PDE problem
%       -div(c*grad(u))+a*u=f where the problem geometry and boundary
%       conditions given by G and B.
%
%       The solution u is represented as the MATLAB column vector U.
%       See ASSEMPDE for details.
%
%       G describes the geometry of the PDE problem. G can
%       either be a Decomposed Geometry Matrix or the name of Geometry
%       MATLAB-file. See either DECSG or PDEGEOM for details.
%
%       B describes the boundary conditions of the PDE problem.  B
%       can either be a Boundary Condition Matrix or the name of Boundary
%       MATLAB-file. See PDEBOUND for details.
%
%       The adapted triangular mesh of the PDE problem is given by the triangle
%       data P, E, and T.  See either INITMESH or PDEGEOM for details.
%
%       The coefficients C, A and F of the PDE problem can
%       be given in a wide variety of ways.  See ASSEMPDE for details.
%
%       Valid property/value pairs include
%
%       Prop.   Value/{Default}         Description
%       ----------------------------------------------------------------------
%       Maxt    Positive integer {Inf}  Maximum number of new triangles
%       Ngen    Positive integer {10}   Maximum number of triangle generations
%       Mesh    P1, E1, T1              Initial mesh
%       Tripick {pdeadworst}|pdeadgsc   Triangle selection method
%       Par     Numeric {0.5}           Function parameter
%       Rmethod {longest}|regular       Triangle refinement method
%       Nonlin  on | off                Use nonlinear solver
%       Toln    numeric {1e-3}          Nonlinear tolerance
%       Init    string|numeric          Nonlinear initial solution value
%       Jac     {fixed}|lumped|full     Nonlinear solver Jacobian calculation
%       Norm    Numeric {Inf}           Nonlinear solver residual norm
%       MesherVersion   {preR2013a}|R2013a      Meshing algorithm used
%
%       Par is passed to the tripick function. Normally it is used as
%       tolerance of how well the solution fits the equation. No more than
%       Ngen successive refinements are attempted. Refinement is also
%       stopped when the number of triangles in the mesh exceeds the
%       Maxt.
%
%       P1, E1, and T1 are the input mesh data. This triangle mesh is used
%       as a starting mesh for the adaptive algorithm. If no initial mesh
%       is provided, the result of a call to INITMESH with no options is
%       used as initial mesh.
%
%       The triangle pick method is a user-definable triangle selection
%       method.  Given the error estimate computed by the function PDEJMPS,
%       the triangle pick method selects the triangles to be refined
%       in the next triangle generation. The function is called using the
%       arguments P, T, CC, AA, FF, U, ERRF, and PAR.  P and T represent the
%       current generation of triangles, CC, AA, FF are the current
%       coefficients for the PDE problem, expanded to triangle midpoints,
%       U is the current solution, ERRF is the computed error estimate,
%       and PAR, the function parameter, given to ADAPTMESH as optional
%       argument. The matrices CC, AA, FF, and ERRF all have NT columns,
%       where NT is the current number of triangles.
%       The number of rows in CC, AA, and FF are exactly the same as the
%       input arguments C, A, and F. ERRF has one row for each equation
%       in the system. There are two standard triangle selection methods
%       in the PDE Toolbox - PDEADWORST and PDEADGSC.
%       PDEADWORST selects triangles where ERRF exceeds a fraction
%       (default: 0.5) of the worst value. PDEADGSC selects triangles
%       using a relative tolerance criterion.
%
%       The refinement method is either 'longest' or 'regular'.
%       See REFINEMESH for details.
%
%       Also nonlinear PDE problems can be solved by the adaptive algorithms.
%       For nonlinear PDE problems, the 'Nonlin' parameter must be set
%       to 'on'. The nonlinear tolerance Toln and nonlinear initial value
%       U0 are passed to the nonlinear solver. See PDENONLIN for details.
%
%       The MesherVersion parameter selects the specific method used to
%       produce a mesh. Option R2013a selects a Delaunay-based algorithm that is
%       faster and more robust at dealing with non-standard geometries. The
%       default option, preR2013a, selects the method used in previous versions
%       of the toolbox.
%
%       See also ASSEMPDE, PDEBOUND, PDEGEOM, INITMESH, REFINEMESH, PDENONLIN

%       Copyright 1994-2012 The MathWorks, Inc.

alfa=0.15;
beta=0.15;
mexp=1;
mesh=0;
nonl=0;
gotu=0;

% Default values
Tripick='pdeadworst'; %#function pdeadworst
Rmethod='longest';
Toln=1e-4;
Ngen=10;
Maxt=Inf;
Par=0.5;
Jac='fixed';
norm=Inf;
mesherVer = 'preR2013a';

k=1;
noptarg=nargin-5;

while k<=noptarg
    Name=eval(['p' int2str(k)]);
    if ~ischar(Name)
        error(message('pde:adaptmesh:ParamNotString'))
    elseif size(Name,1)~=1,
        error(message('pde:adaptmesh:ParamNumRowsOrEmpty'))
    end
    Name=lower(Name);
    if strcmp(Name,'mesh')
        if noptarg-k<3
            error(message('pde:adaptmesh:MeshNumValues'))
        end
        k=k+1;
        p=eval(['p' int2str(k)]);
        k=k+1;
        e=eval(['p' int2str(k)]);
        k=k+1;
        t=eval(['p' int2str(k)]);
        if ischar(p) || ischar(e) || ischar(t)
            error(message('pde:adaptmesh:MeshNotNumeric'));
        end
        mesh=1;
    elseif strcmp(Name,'tripick')
        if noptarg-k<1
            error(message('pde:adaptmesh:TripickNoValue'))
        end
        k=k+1;
        Tripick=eval(['p' int2str(k)]);
        if ~ischar(Tripick)
            error(message('pde:adaptmesh:TripickNotString'))
        end
    elseif strcmp(Name,'rmethod')
        if noptarg-k<1
            error(message('pde:adaptmesh:RmethodNoValue'))
        end
        k=k+1;
        Rmethod=eval(['p' int2str(k)]);
        if ~ischar(Rmethod)
            error(message('pde:adaptmesh:RmethodNotString'))
        end
    elseif strcmp(Name,'toln')
        if noptarg-k<1
            error(message('pde:adaptmesh:TolnNoValue'))
        end
        k=k+1;
        Toln=eval(['p' int2str(k)]);
        if ischar(Toln)
            error(message('pde:adaptmesh:TolnNotNumeric'))
        end
    elseif strcmp(Name,'ngen')
        if noptarg-k<1
            error(message('pde:adaptmesh:NgenNoValue'))
        end
        k=k+1;
        Ngen=eval(['p' int2str(k)]);
        if ischar(Ngen)
            error(message('pde:adaptmesh:NgenNotNumeric'))
        end
    elseif strcmp(Name,'maxt')
        if noptarg-k<1
            error(message('pde:adaptmesh:MaxtNoValue'))
        end
        k=k+1;
        Maxt=eval(['p' int2str(k)]);
        if ischar(Maxt)
            error(message('pde:adaptmesh:MaxtNotNumeric'))
        end
    elseif strcmp(Name,'par')
        if noptarg-k<1
            error(message('pde:adaptmesh:ParNoValue'))
        end
        k=k+1;
        Par=eval(['p' int2str(k)]);
        if ischar(Par)
            error(message('pde:adaptmesh:ParNotNumeric'))
        end
    elseif strcmp(Name,'nonlin')
        if noptarg-k<1
            error(message('pde:adaptmesh:NonlinNoValue'))
        end
        k=k+1;
        Nonlin=eval(['p' int2str(k)]);
        if ~ischar(Nonlin)
            error(message('pde:adaptmesh:NonlinNotString'))
        end
        Nonlin=lower(Nonlin);
        if strcmp(Nonlin,'on')
            nonl=1;
        elseif ~strcmp(Nonlin,'off')
            error(message('pde:adaptmesh:NonlinInvalidString'))
        end
    elseif strcmp(Name,'init')
        if noptarg-k<1
            error(message('pde:adaptmesh:InitNoValue'))
        end
        k=k+1;
        u=eval(['p' int2str(k)]);
        gotu=1;
    elseif strcmp(Name,'jac')
        if noptarg-k<1
            error(message('pde:adaptmesh:JacNoValue'))
        end
        k=k+1;
        Jac=eval(['p' int2str(k)]);
        if ~ischar(Jac)
            error(message('pde:adaptmesh:JacNotString'))
        end
        Jac=lower(deblank(Jac));
        if ~(strcmp(Jac,'fixed') || strcmp(Jac,'lumped') || strcmp(Jac,'full'))
            error(message('pde:adaptmesh:JacInvalidString'))
        end
    elseif strcmp(Name,'norm')
        if noptarg-k<1
            error(message('pde:adaptmesh:NormNoValue'))
        end
        k=k+1;
        norm=eval(['p' int2str(k)]);
        if ischar(norm)
            error(message('pde:adaptmesh:NormNotNumeric'))
        end
    elseif strcmp(Name,'mesherversion')
        if noptarg-k<1
            error(message('pde:adaptmesh:MesherVersionNoValue'))
        end
        k=k+1;
        mesherVer = eval(['p' int2str(k)]);
        if(~ischar(mesherVer))
            error(message('pde:adaptmesh:MesherVersionNotString'))
        end
        if ~(strcmp(mesherVer,'preR2013a') || strcmp(mesherVer,'R2013a'))
            error(message('pde:adaptmesh:MesherVersionInvalidString'))
        end
    else
        error(message('pde:adaptmesh:InvalidOption', Name))
    end
    k=k+1;
end
h = waitbar(0,'1','Name',['Generating adaptMesh ', num2str(Maxt)],...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
if ~mesh
    [p,e,t]=initmesh(g, 'MesherVersion', mesherVer);
end
np=size(p,2);
gen=0;
while 1,
    waitbar(size(t,2)/Maxt,h,sprintf('Number of triangles: %g',size(t,2)))
    %   fprintf('Number of triangles: %g\n',size(t,2))
    if nonl
        if gotu
            u=pdenonlin(b,p,e,t,c,a,f,'jacobian',Jac,'U0',u,'tol',Toln,'norm',norm);
        else
            u=pdenonlin(b,p,e,t,c,a,f,'jacobian',Jac,'tol',Toln,'norm',norm);
        end
        gotu=1;
    else
        u=assempde(b,p,e,t,c,a,f);
    end
    
    if any(isnan(u))
        error(message('pde:adaptmesh:NaNinSolution'))
    end
    
    % Expand values
    [cc,aa,ff]=pdetxpd(p,t,u,c,a,f);
    
    errf=pdejmps(p,t,cc,aa,ff,u,alfa,beta,mexp);
    
    i=feval(Tripick,p,t,cc,aa,ff,u,errf,Par);
    
    if isempty(i),
        waitbar(0,h,sprintf('Adaption completed. %g',size(t,2)))
        fprintf('Adaption completed.\n')
        delete(h) 
        code=1;
        break;
    elseif size(t,2)>Maxt
        waitbar(1,h,sprintf('Maximum number of triangles obtained %g',size(t,2)))
        fprintf('Maximum number of triangles obtained.\n');
        delete(h) 
        code=2;
        break
    elseif gen>=Ngen,
        waitbar(1,h,sprintf('Maximum number of triangles obtained %g',size(t,2)))
        fprintf('\nMaximum number of refinement passes obtained.\n');
        delete(h) 
        code=3;
        break
    elseif getappdata(h,'canceling')
        delete(h)
        code=4;
        break
    end
    
    tl=i';
    % Kludge: tl must be a column vector
    if size(tl,1)==1,
        tl=[tl;tl];
    end
    
    u=reshape(u,np,length(u)/np);
    [p,e,t,u]=refinemesh(g,p,e,t,u,tl,Rmethod);
    u=u(:);
    np=size(p,2);
    gen=gen+1;
end


