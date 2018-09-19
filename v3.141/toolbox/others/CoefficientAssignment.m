classdef (Sealed) CoefficientAssignment  < pde.EquationAssignment & matlab.mixin.internal.CompactDisplay 
% CoefficientAssignment Specify all PDE coefficients over a domain or subdomain
%     The PDE toolbox can solve equations of the form:
%  
%              m*d^2u/dt^2 + d*du/dt - div(c*grad(u)) + a*u = f
%  
%     and the corresponding eigenvalue equation of the form:
%  
%                   - div(c*grad(u)) + a*u = lamda*d*u
%     or
%                   - div(c*grad(u)) + a*u = (lamda^2)*m*u
%  
%     The equation to solve is defined in terms of the coefficients m, d, c, 
%     a, f. This method creates an object representing the coefficients 
%     in a domain or subdomain and appends the object to the EquationCoefficients 
%     property.
%  
%     You can define more than one set of equation coefficients. For example, 
%     for a 2-D geometry that has two subdomains and equation coefficients
%     that differ in each subdomain. However, the same coefficient terms must
%     be present in each equation - they just differ in value.
%     Note: Subdomains are not currently supported in 3-D.
%
%     Instances of this class can only be created by calling the 
%     specifyCoefficients method of the PDEModel class. 
%
% See also pde.PDEModel, pde.PDEModel/specifyCoefficients

% Copyright 2015 The MathWorks, Inc.


properties (SetAccess = private)
% RegionType - Type of geometric region the coefficients are assigned to.
%    A string specifying the type of geometric domain the coefficients
%    are assigned to. This string has two possible values: 'cell'
%    or 'face'. 
    RegionType;
end

properties
% RegionID - ID of the geometric regions the coefficients are assigned to.
%    For 2-D each RegionID satisfies 0 < RegionID(j) < NumFaces in the 
%    geometry. For 3-D each RegionID satisfies 0 < RegionID(j) < NumCells 
%    in the geometry.  
    RegionID;

% m - Second-order derivative coefficient
%     A numeric scalar, vector or matrix, or function
%     handle representing the m coefficient of the PDE.
     m;

% d - First-order derivative coefficient
%     A numeric scalar, vector or matrix, or function
%     handle representing the d coefficient of the PDE.
     d;

% c - Grad term coefficient        
%     A numeric scalar, vector or matrix, or function
%     handle representing the c coefficient of the PDE.
     c;

% a - Dependent variable u coefficient   
%     A numeric scalar, vector or matrix, or function
%     handle representing the a coefficient of the PDE.
     a;

% f - Right-hand side source
%     A numeric scalar, vector or matrix, function
%     handle representing the source term of the PDE.
     f;
end

methods(Hidden=true, Access={?pde.PDEModel})
    function obj=CoefficientAssignment(car, varargin) 
      obj.RecordOwner = car; 
      parser = inputParser;
      parser.PartialMatching=false; % Clash between 'face' and 'f'     
      parser.addParameter('face', []);     
      parser.addParameter('cell', []);     
      parser.addParameter('m', []);                 
      parser.addParameter('d', []);
      parser.addParameter('c', []);
      parser.addParameter('a', []);    
      parser.addParameter('f', []); 
      parser.addParameter('SystemSize', 1); 
      parser.parse(varargin{:});
           
      numdims = 2;
      if ~isempty(parser.Results.face)
              obj.RegionType = 'face';
              obj.RegionID = parser.Results.face;
      elseif ~isempty(parser.Results.cell)
             obj.RegionType = 'cell';
             obj.RegionID = parser.Results.cell; 
             numdims = 3;
      end       
      systemsize = parser.Results.SystemSize;

      obj.m = parser.Results.m;
      obj.d = parser.Results.d;
      obj.c = parser.Results.c;   
      obj.a = parser.Results.a;
      obj.f = parser.Results.f;            
      obj.checkAllMatrixCoefSizes(systemsize, numdims);  
      obj.checkFcnHdlArgCounts(systemsize, numdims);
      obj.checkMandD();
      obj.checkSparseness();
      % obj.checkNumericComplexity();
    end
end

methods
    
     function set.RegionType(self, rtype)     
        self.ValidateRegionType(rtype);     
        self.RegionType = rtype;       
    end
    
    function set.RegionID(self, rids)      
        self.ValidateRegionID(rids);      
        self.RegionID = rids;        
    end        
    
    function set.m(self, coef)    
      self.CoefPrecheck(coef);     
      self.m = coef;     
    end
    
    function set.d(self, coef) 
      self.CoefPrecheck(coef);
      self.d = coef;   
    end
    
    function set.c(self, coef)  
      self.CoefPrecheck(coef);      
      self.c = coef;   
    end
    
    function set.a(self, coef)    
      self.CoefPrecheck(coef);  
      self.a = coef;     
    end
    
    function set.f(self, coef)  
      self.CoefPrecheck(coef);
      self.f = coef;  
    end

end

methods(Hidden=true, Access = {?pde.CoefficientAssignmentRecords})
    function tf=sameCoefficients(self,other,varargin)
      if ~coefficientsMatch(self, other)
          tf = false;
          return   
      end
      if isempty(varargin)
             tf = (isequal(self.m, other.m) && ...
                   isequal(self.d, other.d) && ...
                   isequal(self.c, other.c) && ...
                   isequal(self.a, other.a) && ...
                   isequal(self.f, other.f));                      
      else
          cname = varargin{1};
          switch cname
            case 'm'
                tf = isequal(self.m, other.m);
            case 'd'
                tf = isequal(self.d, other.d);
            case 'c'
                tf = isequal(self.c, other.c); 
            case 'a'
                tf = isequal(self.a, other.a);    
            case 'f'
                tf = isequal(self.f, other.f);                                  
          end
      end
    end   
    
    function tf=numericCoefficients(self,varargin) 
      if isempty(varargin)
          
             tf = ~any(isa(self.m, 'function_handle')  || ...
                      isa(self.d, 'function_handle')  || ...
                      isa(self.c, 'function_handle')  || ...
                      isa(self.a, 'function_handle')  || ...
                      isa(self.f, 'function_handle'));                               
      else
          cname = varargin{1};
          switch cname
            case 'm'
                tf = ~isa(self.m, 'function_handle');
            case 'd'
                tf = ~isa(self.d, 'function_handle');
            case 'c'
                tf = ~isa(self.c, 'function_handle');
            case 'a'
                tf = ~isa(self.a, 'function_handle');   
            case 'f'
                tf = ~isa(self.f, 'function_handle');                            
          end
      end
    end 
            
    
    function tf = coefficientsMatch(self, other)
        tf = true;
        tf = tf & (self.mDefined() == other.mDefined());       
        tf = tf & (self.dDefined() == other.dDefined());         
        tf = tf & (self.cDefined() == other.cDefined());             
%         tf = tf & (self.aDefined() == other.aDefined());                    
%         tf = tf & (self.fDefined() == other.fDefined());           
    end
    
    function performSolverPrecheck(self, systemsize, numfaces, numcells)
         ndims = 2;       
         if strcmp(self.RegionType, 'face')
            if any(self.RegionID > numfaces)
                error(message('pde:pdeCoefficientSpecification:invalidFaceIndexPresolve'));
            end
         else
            if any(self.RegionID > numcells)
                error(message('pde:pdeCoefficientSpecification:invalidCellIndexPresolve'));
            end 
            ndims = 3;
         end  
         checkAllMatrixCoefSizes(self, systemsize, ndims); 
         checkMandD(self);
         checkSparseness(self);
         checkFcnHdlArgCounts(self, systemsize, ndims)   
         % checkNumericComplexity(self); 
     end                
    
end


methods(Hidden=true, Access = private)        
    % checkAllMatrixCoefSizes - check the size of all coefficients that
    % are defined by a matrix.
    function checkAllMatrixCoefSizes(self, systemsize, ndims)
       self.checkMCoefSize(self.m, systemsize);
       if ~(self.mDefined() && self.dDefined())
          self.checkDCoefSize(self.d, systemsize);
       end
       self.checkCCoefSize(self.c, systemsize,ndims);
       self.checkACoefSize(self.a, systemsize);
       self.checkFCoefSize(self.f, systemsize);        
    end
    function checkSparseness(self)
        sparsecoef = issparse(self.m) | issparse(self.c) | issparse(self.a) | issparse(self.f);
        if ~(self.mDefined() && self.dDefined())
           sparsecoef = sparsecoef |  issparse(self.d);
        end
        if sparsecoef
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueSparse'));   
        end        
    end  
    
    function checkFcnHdlArgCounts(self, systemsize, ndims)        
        if self.mDefined() 
            self.checkCoefFcnHdlArgCounts(self.m, systemsize, ndims);              
        end                       
        if self.dDefined() 
            self.checkCoefFcnHdlArgCounts(self.d, systemsize, ndims);         
        end               
        if self.cDefined()
            self.checkCoefFcnHdlArgCounts(self.c, systemsize, ndims);     
        end               
        if self.aDefined() 
            self.checkCoefFcnHdlArgCounts(self.a, systemsize, ndims);      
        end               
        if self.fDefined() 
            self.checkCoefFcnHdlArgCounts(self.f, systemsize, ndims);    
        end             
    end
    
end

methods(Hidden=true, Access = {?pde.CoefficientAssignmentRecords})
     function tf = mDefined(self)
        tf = self.coefDefined(self.m); 
     end
    
     function tf = dDefined(self)
        tf = self.coefDefined(self.d); 
     end

     function tf = cDefined(self)
        tf = self.coefDefined(self.c); 
     end

     function tf = aDefined(self)
        tf = self.coefDefined(self.a); 
     end

     function tf = fDefined(self)
        tf = self.coefDefined(self.f); 
     end
        
    function checkMandD(self)        
        if (self.mDefined() && self.dDefined())
           if isa(self.d,'function_handle') || isscalar(self.d) || isvector(self.d)
               error(message('pde:pdeCoefficientSpecification:dmustBeMatrix'));  
           end
        end       
    end    
%     function checkNumericComplexity(self) 
%         import pde.CoefficientAssignment.*
%         realcoef = true(1,0);  
%         if self.mDefined() 
%             realcoef(end+1) = coefIsComplexNumeric(self.m);              
%         end                       
%         if (self.dDefined() && isnumeric(self.d))
%             realcoef(end+1) = coefIsComplexNumeric(self.d);         
%         end               
%         if (self.cDefined() && isnumeric(self.c))
%             realcoef(end+1) = coefIsComplexNumeric(self.c);     
%         end               
%         if (self.aDefined() && isnumeric(self.a))
%             realcoef(end+1) = coefIsComplexNumeric(self.a);      
%         end               
%         if (self.fDefined() && isnumeric(self.f))
%             realcoef(end+1) = coefIsComplexNumeric(self.f);    
%         end       
%         if all(realcoef == false) || all(realcoef == true)
%             return
%         else
%             error(message('pde:pdeCoefficientSpecification:complexRealMix'));  
%         end               
%     end    
    
   function tf = hasComplexCoefficient(self, loc, state)  
        import pde.CoefficientAssignment.*
        tf = false(5,1);        
        tf(1) = coefIsComplexNumericOrFcnHdl(self.m, loc, state);
        tf(2) = coefIsComplexNumericOrFcnHdl(self.d, loc, state);
        tf(3) = coefIsComplexNumericOrFcnHdl(self.c, loc, state);
        tf(4) = coefIsComplexNumericOrFcnHdl(self.a, loc, state);
        tf(5) = coefIsComplexNumericOrFcnHdl(self.f, loc, state);                      
   end    
   
end


methods(Static, Hidden=true)
    function preDelete(self,~)
        if isvalid(self.RecordOwner)
            self.RecordOwner.delistCoefficientAssignment(self);
        end
    end    
end

methods(Static, Access = private)
    function tf = coefDefined(coef)
       tf = (isnumeric(coef) && ~(isscalar(coef) && coef == 0) || isa(coef,'function_handle'));
    end
    function tf = coefIsComplexNumeric(coef)
         tf = false;
         if ~pde.CoefficientAssignment.coefDefined(coef)
             return 
         end
         if isnumeric(coef)
            if ~isreal(coef)
                tf = true;          
            end        
         end              
    end         
    function tf = coefIsComplexNumericOrFcnHdl(coef, loc, state)
         tf = false;
         if ~pde.CoefficientAssignment.coefDefined(coef)
             return 
         end
         if isnumeric(coef)
            if ~isreal(coef)
                tf = true;          
            end 
         else % isa(coef, 'function_handle')
             res = coef(loc, state);
             if ~isreal(res)
                tf = true;          
             end 
         end              
    end        
    function ok=ValidateRegionType(rgntype)   
      nc = numel(rgntype);
      if ~(strncmpi(rgntype,'face',nc) || strncmpi(rgntype,'cell',nc))
        error(message('pde:pdeCoefficientSpecification:invalidRegionType'));   
      end
      ok = true;      
    end
    function ok=ValidateRegionID(rgnid)
      % Must be real(non-complex), full, natural number.        
      if ~isreal(rgnid) || ~all(rgnid(:) > 0) || issparse(rgnid) || any(mod(rgnid(:),1)~=0)
        error(message('pde:pdeCoefficientSpecification:invalidRegionID'));   
      end      
      ok = true; 
    end
    function ok=CoefPrecheck(coef)
      ok =false;
      if isfloat(coef)                   
        if any(isnan(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueNaN'));
        elseif any(isinf(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueInf'));
        elseif(isempty(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueEmpty'));
        end
      elseif ~isa(coef, 'function_handle')            
          if pde.CoefficientAssignment.stringCoefMissedDef(coef)
            error(message('pde:pdeCoefficientSpecification:missedCoefValue'));   
          elseif ischar(coef)
             error(message('pde:pdeCoefficientSpecification:invalidCoefValueString')); 
          else
            error(message('pde:pdeCoefficientSpecification:invalidCoefValue'));       
          end                               
      end      
      ok = true; 
    end
    
    function checkMCoefSize(mcoef, systemsize)
          if isscalar(mcoef) && ~isa(mcoef, 'function_handle') && mcoef == 0
            return
          elseif isa(mcoef, 'function_handle')
            return  
          end
         mveclen = numel(mcoef);
         lengthok = (mveclen == 0 || mveclen == 1 || mveclen == systemsize || ...
                     mveclen == (systemsize*(systemsize+1)/2) || mveclen == systemsize*systemsize);          
         if ~(isvector(mcoef)  && iscolumn(mcoef) && lengthok)      
            error(message('pde:pdeCoefficientSpecification:invalidAMatrixSize'));   
        end
    end
    function checkDCoefSize(dcoef, systemsize)
         if isscalar(dcoef) && ~isa(dcoef, 'function_handle') && dcoef == 0
            return
         elseif isa(dcoef, 'function_handle')
             return
         end
         dveclen = numel(dcoef);
         lengthok = (dveclen == 0 || dveclen == 1 || dveclen == systemsize || ...
                     dveclen == (systemsize*(systemsize+1)/2) || dveclen == systemsize*systemsize);          
         if ~(isvector(dcoef)  && iscolumn(dcoef) && lengthok)      
            error(message('pde:pdeCoefficientSpecification:invalidDMatrixSize'));   
         end
    end
    function checkCCoefSize(ccoef, systemsize, ndims)
        if isscalar(ccoef) && ~isa(ccoef, 'function_handle') && ccoef == 0
            return
        elseif isa(ccoef, 'function_handle')
            return
        end
        cveclen = numel(ccoef);      
        if ndims == 2
             lengthok = (cveclen == 1 || cveclen == 2 || cveclen == 3 || ...
                         cveclen == 4 || cveclen == systemsize || ...
                         cveclen == 2*systemsize || cveclen == 3*systemsize || ...
                         cveclen == 4*systemsize || cveclen == 2*systemsize*((2*systemsize)+1)/2 || ...
                         cveclen == 4*(systemsize^2));
             if ~(isvector(ccoef)  && iscolumn(ccoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end                    
        else
            lengthok = (cveclen == 1 || cveclen == 3 || cveclen == 6 || ...
                         cveclen == 9 || cveclen == systemsize || ...
                         cveclen == 3*systemsize || cveclen == 6*systemsize || ...
                         cveclen == 9*systemsize || cveclen == 3*systemsize*((3*systemsize)+1)/2 || ...
                         cveclen == 9*(systemsize^2));
             if ~(isvector(ccoef)  && iscolumn(ccoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end               
        end                      
    end
    function checkACoefSize(acoef, systemsize) 
         if isscalar(acoef) && ~isa(acoef, 'function_handle') && acoef == 0
            return
         elseif isa(acoef, 'function_handle') 
            return
         end
         aveclen = numel(acoef);        
         lengthok = (aveclen == 0 || aveclen == 1 || aveclen == systemsize || ...
                     aveclen == (systemsize*(systemsize+1)/2) || aveclen == systemsize*systemsize);          
         if ~(isvector(acoef)  && iscolumn(acoef) && lengthok)      
            error(message('pde:pdeCoefficientSpecification:invalidAMatrixSize'));   
        end
    end   
    function checkFCoefSize(fcoef, systemsize)
         if isscalar(fcoef) && ~isa(fcoef, 'function_handle') && fcoef == 0
            return
         elseif isa(fcoef, 'function_handle')
            return 
         end
        if ~(isvector(fcoef)  && iscolumn(fcoef) && numel(fcoef) == systemsize)      
            error(message('pde:pdeCoefficientSpecification:invalidFMatrixSize'));   
        end
    end
    
    %
    % stringCoefMissedDef
    % Returns true if the coefficient specification is a char of length 1
    % and is either 'm', 'd', 'c','a', 'f'.
    % This detects a typo that's not uncommon;  e.g. 
    % (..., 'm', 1, 'd', 2, 'c', 'a', 'f', 0)
    function tf = stringCoefMissedDef(coef)
        tf = false;
        if ischar(coef) && numel(coef) == 1
            tf = strcmp('m',coef) || strcmp('d',coef) || strcmp('c',coef) ...
                 || strcmp('a',coef) || strcmp('f',coef);
        end
    end           
         
    
%
%   Coefficient function-handle argument checks
%   
    function checkCoefFcnHdlArgCounts(coef, systemsize, ndims)
        if ~isa(coef, 'function_handle')
            return
        end
        location.x = 0;
        location.y = 0;          
        location.z = 0;
        location.subdomain=1;  
        state.u = zeros(systemsize, 1);
        state.ux = zeros(systemsize, 1);
        state.uy = zeros(systemsize, 1);
        state.uz = zeros(systemsize, 1);
        state.time = 0;  
        try 
            coef(location, state);
        catch
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));   
        end
        
        try 
            res = coef(location, state);
        catch
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));   
        end
                     
        if ndims == 2
           p = [0 1 0; 0 0 1];
           t = [1; 2; 3];          
        else            
           p = [0 1 0 0; 0 0 1 0; 0 0 0 1]; 
           t = [1; 2; 3; 4];           
        end
        u = 0;
        time = 0;
        threwerr = false; 
        try 
            res = coef(p,t,u,time);
        catch
            threwerr = true;
        end        
        if ~threwerr
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));    
        end        
    end  
    
end    
   properties (Hidden = true, Access='private')
      RecordOwner;
   end   
end
