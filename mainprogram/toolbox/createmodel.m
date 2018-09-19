% create model, no mesh
function varargout=createmodel(considerheating,varargin)
global modc 
if length(modc.rhe)==2 %crust and mantle two layer
    %% interseismic model
    g=@crustmantlethreeface;
    model = createpde;
    geometryFromEdges(model,g);
    % Boundary Conditions
    applyBoundaryCondition(model,'edge',[1,11],'u',modc.v0);% constant velocity boundary
    applyBoundaryCondition(model,'edge',[7,8,9],'u',0);
    % boundary condition matrix
    %% coseismic model
    %     gn=@crustmantletwoface;
    gn=@fivepointrect;
    modeln = createpde;
    geometryFromEdges(modeln,gn);
    % Boundary Conditions
    applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','on');
%     applyBoundaryCondition(modeln,'edge',2,'u',0); % lock fault, not lock far field
    applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
%     applyBoundaryCondition(modeln,'edge',[2,4,6,8],'u',0);
    varargout{1}=model;
    varargout{2}=g;
    varargout{3}=modeln;
    varargout{4}=gn;
    if considerheating
        considerheating=true;
        %% temperature model
        gt=@crustmantlethreeface;
        modelt = createpde;
        geometryFromEdges(modelt,gt);
        % Boundary Conditions
        D=varargin{1};
        thermalgradient = @(region,state) geothermal(region.y,D);
        applyBoundaryCondition(modelt,'edge',[5,6],'u',modc.mod.Ts); % surface temperature
        mther=modc.ther{2};
        cther=modc.ther{1};
        applyBoundaryCondition(modelt,'edge',10,'u',modc.mod.Tb); % constant base temperature
        applyBoundaryCondition(modelt,'edge',[1,11],'u',thermalgradient); % far field temperature
        
        specifyCoefficients(modelt,'m',0,'d',cther.rho*cther.cp,'c',cther.k,'a',0,'f',0,'Face',[1,2]); % crust
        specifyCoefficients(modelt,'m',0,'d',mther.rho*mther.cp,'c',mther.k,'a',0,'f',0,'Face',3); % mantle
        varargout{5}=modelt;
        varargout{6}=gt;
    end
elseif length(modc.rhe)==1
    %% interseismic model
    g=@rectwithsubdomain;
    model = createpde;
    geometryFromEdges(model,g);
    % Boundary Conditions
    applyBoundaryCondition(model,'edge',1,'u',modc.v0);% constant velocity boundary
    applyBoundaryCondition(model,'edge',[7,8],'u',0);
    % boundary condition matrix
    %% coseismic model
    gn=@fivepointrect;
    modeln = createpde;
    geometryFromEdges(modeln,gn);
    % Boundary Conditions
    applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','off');
    applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
    varargout{1}=model;
    varargout{2}=g;
    varargout{3}=modeln;
    varargout{4}=gn;
    if considerheating
        considerheating=true;
        %% temperature model
        gt=@rectwithsubdomain;
        modelt = createpde;
        geometryFromEdges(modelt,gt);
        % Boundary Conditions
        thermalgradient = @(region,state)0.025*region.y+modc.mod.Ts;
        applyBoundaryCondition(modelt,'edge',[5,6],'u',modc.mod.Ts); % surface temperature
        applyBoundaryCondition(modelt,'edge',2,'g',modc.mod.dT*ther.k); % constant moho heat flow
        applyBoundaryCondition(modelt,'edge',1,'u',thermalgradient); % far field temperature
        specifyCoefficients(modelt,'m',0,'d',rho*ther.cp,'c',ther.k,'a',0,'f',0);
        varargout{5}=modelt;
        varargout{6}=gt;
    end
end
%% apply mesh to mod
specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0);
model.Mesh=node2mesh(modc.Nodes,g);
% jigglemesh
[p,e,t]=meshToPet(model.Mesh);
p = jigglemesh(p,e,t,'opt','mean','iter',inf); 
modc.Nodes=p;
model.Mesh=node2mesh(p,g);
modeln.Mesh=node2mesh(p,gn);
if considerheating
    modelt.Mesh=node2mesh(p,gt);
end
end