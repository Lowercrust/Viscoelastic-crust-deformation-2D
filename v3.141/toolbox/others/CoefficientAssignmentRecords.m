classdef (Sealed) CoefficientAssignmentRecords  < handle & matlab.mixin.internal.CompactDisplay
% CoefficientAssignmentRecords  Assignments of the equation coefficients
%    The CoefficientAssignmentRecords holds a record of all the assignments of 
%    equation coefficients across the geometric domain. This may be a single
%    assignment that represents the same field equations throughout the domain
%    or multiple assignments where the field equations are different across
%    various subdomains. This object cannot be created directly, it can only 
%    be created by the PDEModel class.
%
% CoefficientAssignmentRecords methods:
%    globalCoefficients - True if consistent across the geometric domain
%    findCoefficients   - Find coefficients assignment for a geometric region
%
% CoefficientAssignmentRecords properties:
%    CoefficientAssignments - Vector of coefficient assignments
%
% See also pde.PDEModel, pde.PDEModel/specifyCoefficients

% Copyright 2015 The MathWorks, Inc.

properties
    
% CoefficientAssignments - Vector of coefficient assignments
%    A vector containing the coefficient assignments to the geometric domain
%    or subdomains. Each entry in the vector is a pde.CoefficientAssignment object. 
%    This object defines the coefficients of the PDE that should be imposed
%    on a designated domain or subdomain. Coefficients and host domains are 
%    defined using the speficyCoefficients method of the PDEModel.
    CoefficientAssignments;

end

methods
   % Method declaration
   tf = globalCoefficients(self, varargin)    
   ca = findCoefficients(self, varargin)   
   tf = coefficientsSpecified(self, varargin)  
end

methods(Hidden=true, Access={?pde.PDEModel})
    function obj=CoefficientAssignmentRecords(pdem) 
       obj.CoefficientAssignments = [];
       obj.ParentPdemodel = pdem;
    end  
    
    function tf = allCoefficientsNumeric(self)
        tf = true;
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
        end       
        nf = self.ParentPdemodel.Geometry.NumFaces;      
        if ~self.ParentPdemodel.IsTwoD
            ca = self.CoefficientAssignments(end);
            tf = numericCoefficients(ca);
        else 
            for i = 1:nf
                ca = self.findCoefficients('face',i);
                if ~numericCoefficients(ca)
                    tf = false;
                    return
                end
            end
        end        
    end
    
    % Check that each subdomain has a coefficient assignment
    % also perform a basic sanity check on the coefficient
    % in the event that it was edited since assigned.
    % If both m and d are defined, then d must be a matrix.
    function performSolverPrechecks(self)
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
        end
        if ~consistentCoefficientDefinition(self)
            error(message('pde:pdeCoefficientSpecification:inconsistentCoefPresence')); 
        end
        
%        [location, state] = getSampleLocState(self)        
%         if ~consistentComplexityDefinition(self, location, state)
%             error(message('pde:pdeCoefficientSpecification:complexRealMix'));  
%         end
        
        syssize = self.ParentPdemodel.PDESystemSize;
        nf = self.ParentPdemodel.Geometry.NumFaces;
        nc = 0;
        if ~self.ParentPdemodel.IsTwoD
            nc = self.ParentPdemodel.Geometry.NumCells;
            ca = self.CoefficientAssignments(end);
            ca.performSolverPrecheck(syssize, nf, nc);   
        else           
            for i = 1:nf
                ca = self.findCoefficients('face',i);
                ca.performSolverPrecheck(syssize, nf, nc);   
            end
        end                         
    end
    
    function coefstruct = packCoefficients(self)  
        [loc, state] = getSampleLocState(self);
        msh = self.ParentPdemodel.Mesh;
        ma = msh.MeshAssociation;
        eidtofid = ma.RegionAssociativity; % Element ID to face ID.
        numelems = size(msh.Elements,2);        
        if self.globalCoefficients()           
            coefstruct.ElementsInSubdomain = (1:numelems)'; 
            coefstruct.NumElementsInSubdomain = numelems;
            ca = self.CoefficientAssignments(end);                     
            thisstruct.m{1} = ca.m;
            thisstruct.d{1} = ca.d;
            thisstruct.c{1} = ca.c;
            thisstruct.a{1} = ca.a;
            thisstruct.f{1} = ca.f;
            coefstruct.Coefficients = thisstruct;            
            coefstruct.IsComplex = hasComplexCoefficient(ca, loc, state);  
        else
            nf = self.ParentPdemodel.Geometry.NumFaces;
            coefstruct.ElementsInSubdomain = zeros(numelems,1)';
            coefstruct.NumElementsInSubdomain = zeros(nf,1); 
            complexcoef = false(5,1);
            endidx = 0;
            for i = 1:nf             
                felems = find(eidtofid==i);
                numfelems = numel(felems);
                startidx = endidx;
                endidx = startidx+numfelems;
                startidx = startidx+1;
                coefstruct.NumElementsInSubdomain(i) = numfelems;
                coefstruct.ElementsInSubdomain(startidx:endidx)=felems;              
                ca = self.findCoefficients('face',i);                                             
                thisstruct.m{i} = ca.m;
                thisstruct.d{i} = ca.d;
                thisstruct.c{i} = ca.c;
                thisstruct.a{i} = ca.a;
                thisstruct.f{i} = ca.f;
                coefstruct.Coefficients = thisstruct;                                               
                complexcoef = complexcoef | hasComplexCoefficient(ca, loc, state);  
            end  
            coefstruct.IsComplex = complexcoef;
        end        
    end    
end    
    
methods(Hidden=true, Access={?pde.PDEModel,?pde.InitialConditionsRecords})    
    function tf = mDefined(self)
       tf = false;
       numassigns = numel(self.CoefficientAssignments);
       if numassigns
         ca = self.CoefficientAssignments(end);
         tf = ca.mDefined();
       end   
    end
    
   function tf = dDefined(self)
       tf = false;
       numassigns = numel(self.CoefficientAssignments);
       if numassigns
         ca = self.CoefficientAssignments(end);
         tf = ca.dDefined();
       end   
   end
    
   function tf = cDefined(self)
       tf = false;
       numassigns = numel(self.CoefficientAssignments);
       if numassigns
         ca = self.CoefficientAssignments(end);
         tf = ca.cDefined();
       end   
    end
    
    function tf = aDefined(self)
       tf = false;
       numassigns = numel(self.CoefficientAssignments);
       if numassigns
         ca = self.CoefficientAssignments(end);
         tf = ca.aDefined();
       end   
    end
    
    function tf = fDefined(self)
       tf = false;
       numassigns = numel(self.CoefficientAssignments);
       if numassigns
         ca = self.CoefficientAssignments(end);
         tf = ca.fDefined();
       end   
    end 
    
    function tf = timeDependent(self)
        tf = (mDefined(self) || dDefined(self));
    end
end
  
methods(Hidden=true, Access=private)
    % Returns true if each coef has consistent definition across all
    % subdomains. For example, m coefficient defined in one region is
    % also defined in all other regions.
    function tf = consistentCoefficientDefinition(self)
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
        end
        tf = globalCoefficients(self);
        if tf
            return
        end
        tf = true;
        nf = self.ParentPdemodel.Geometry.NumFaces;             
        thiscoef = self.findCoefficients('face', 1);
        for i = 2:nf
            othercoef = self.findCoefficients('face', i);           
            if ~thiscoef.coefficientsMatch(othercoef)
                tf = false;
                return;
            end            
        end                   
    end   
    
%     function tf = consistentComplexityDefinition(self, loc, state)
%         tf = true;
%         if ~self.coefficientsSpecified()
%             error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
%         end
%         if globalCoefficients(self);
%             return
%         end      
%         nf = self.ParentPdemodel.Geometry.NumFaces;
%         cxcoef = false(nf,1);     
%         for i = 1:nf
%             thiscoef = self.findCoefficients('face', i);           
%             cxcoef(i) = thiscoef.hasComplexCoefficient(loc, state);                
%         end       
%         tf = (all(cxcoef == true) || all(cxcoef == false));
%     end   
    
    
    function tf = hasComplexCoefficients(self)
        
        [location, state] = getSampleLocState(self);                       
        tf = false;
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
        end     
        if globalCoefficients(self)
            % Just need to check one.
            thiscoef = self.CoefficientAssignments(end);
            tf = any(thiscoef.hasComplexCoefficient(location, state));
            return
        end   
        nf = self.ParentPdemodel.Geometry.NumFaces;                    
        for i = 1:nf
            thiscoef = self.findCoefficients('face', i);           
            if any(thiscoef.hasComplexCoefficient(location, state))
                tf = true;
                return;
            end            
        end                   
    end 
    
    function [location, state] = getSampleLocState(self)
        msh = self.ParentPdemodel.Mesh;
        nodecoords =msh.Nodes(:,1);
        location.x = nodecoords(1);
        location.y = nodecoords(2);
        if numel(nodecoords) == 3
            location.z = nodecoords(3);
        else
            location.z = 0;
        end
        location.subdomain=1;         
        N = self.ParentPdemodel.PDESystemSize;
        state.u(1:N,1)=0;   % state.u(1:NSystemSize, NumLocations)
        state.ux(1:N,1)=0;
        state.uy(1:N,1)=0;
        state.uz(1:N,1)=0;
        state.time=0;        
    end
    
end

methods(Hidden=true, Access={?pde.CoefficientAssignment})
    function delistCoefficientAssignment(self, coeftodelete)
        numcoef = numel(self.CoefficientAssignments);
        for i = 1:numcoef
            thiscoef = self.CoefficientAssignments(i);
            if thiscoef == coeftodelete
                self.CoefficientAssignments(i) = [];
                break
            end
        end  
        numcoef = numel(self.CoefficientAssignments);
        if numcoef == 0 && isvalid(self.ParentPdemodel)
           self.ParentPdemodel.delistCoefficientAssignments(); 
        end
    end 
end

methods(Static, Hidden=true)
    function preDelete(self,~)
        if isvalid(self.ParentPdemodel)
            self.ParentPdemodel.delistCoefficientAssignments();
        end
    end    
end

properties (Hidden = true, SetAccess='private')
    ParentPdemodel;
end  
  
end
