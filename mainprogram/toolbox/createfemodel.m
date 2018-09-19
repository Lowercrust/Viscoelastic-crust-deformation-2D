function femodel=createfemodel(model,method,varargin)
global gpunum
% K A F Q G H R M
% coefstruct = packCoefficientssa(model.EquationCoefficients);
[p,~,t]=meshToPet(model.Mesh);
c=model.EquationCoefficients.CoefficientAssignments.c;
femodel=struct('K',{},'K0',{},'A',{},'F',{},'Q',{},'G',{},'H',{},'R',{},'M',{},'B',{},...
    'Mass',{},'Nu',{},'Or',{},'ud',{},'m',{},'L0',{},'U0',{},'P0',{},'Q0',{},'R0',{});
if ~isempty(varargin)
    fmatrix=varargin{1};
else
    fmatrix=zeros(1,length(t));
end
% if model.IsTimeDependent
    if strcmp(method,'new')
        [femodel(1).Q,femodel(1).G,femodel(1).H,femodel(1).R]=assembleBoundary(model);
        if gpunum~=0
            [femodel(1).K,femodel(1).F]=createGlobalKFGPU(p,t,'K',c,'F',fmatrix);
        else
            [femodel(1).K,femodel(1).F]=createGlobalKFs(p,t,'K',c,'F',fmatrix);
        end
        if model.IsTimeDependent
            coefstruct = packCoefficientssa(model.EquationCoefficients);
            % Mass = formGlobalM2D(self.thePde, self.p, self.t, self.coefstruct,u,time,'m');
            femodel(1).Mass = formGlobalM2D(model, p, t, coefstruct,[],eps,'d'); % independent with time and temperature field
            femodel(1).A = formGlobalM2D(model, p, t, coefstruct,[],0,'a'); % 3.5s
        else
            femodel(1).A=sparse(size(p,2),size(p,2),0);
            femodel(1).Mass=sparse(size(p,2),size(p,2),0);
        end
        [KK, FF, B, ud, MM] = applyBCNullSpacesa(femodel.K,femodel.A,femodel.F,femodel.Q,femodel.G,femodel.H,femodel.R,femodel.Mass);
        femodel.K0=femodel.K;
        femodel.K=KK;
        femodel.F=FF;
        femodel.M=MM;
        femodel.B=B;
        [femodel.Nu,femodel.Or]=pdenullorth(femodel.H);
        femodel.ud=ud;
        femodel.m=femodel.B'*femodel.Mass*femodel.B;
        [L,U,P,Q,R]=lu(femodel.m);
        femodel.L0=L;
        femodel.U0=U;
        femodel.P0=P;
        femodel.Q0=Q;
        femodel.R0=R;
%         G = [1; 3/2; 11/6; 25/12; 137/60];
%         alpha = [-37/200; -1/9; -0.0823; -0.0415; 0];
%         invGa = 1 ./ (G .* (1 - alpha));
%         
%         hinvGak = h * invGa(1);
%         Miter = Mt + hinvGak * femodel.K;
%         [L,U,P,Q,R]=lu(femodel.m);
%         femodel.L1=L;
%         femodel.U1=U;
%         femodel.P1=P;
%         femodel.Q1=Q;
%         femodel.R1=R;
        %full(femodel.Or*((femodel.H*femodel.Or)\femodel.R));
        %     femodel=pde.DynamicDiscretizedPDEModel(self,p,e,t,coefstruct,u0,tlist,tsecondOrder);
    elseif strcmp(method,'update') % update source term
        femodel=varargin{2};
        fmatrix=varargin{1};
        if gpunum~=0
            femodel(1).F=createGlobalKFGPU(p,t,'F',fmatrix);
        elseif gpunum==0
            femodel(1).F=createGlobalKFs(p,t,'F',fmatrix);
        end
%         Fc=null'*((F+G)-Kc*ud);
        femodel(1).F=femodel.Nu'*((femodel(1).F+femodel(1).G)-(femodel.K0+femodel.A+femodel.Q)*femodel.ud);
    end
% elseif ~model.IsTimedepent
%     if strcmp(method,'new')
%         % bmatrix=struct('Q',{},'G',{},'H',{},'R',{});
%         [femodel(1).Q,femodel(1).G,femodel(1).H,femodel(1).R]=assembleBoundary(model);
%         if gpunum~=0 % if there is a gpu device,
%             [femodel(1).K,femodel(1).F]=createGlobalKFGPU(p,t,'K','F',fmatrix); %stiffness matrix
%         else
%             [femodel(1).K,femodel(1).F]=createGlobalKFs(p,t,'K',c,'F',fmatrix); %stiffness matrix
%         end
%         femodel(1).A=sparse(size(p,2),size(p,2),0);
%         femodel(1).M=sparse(size(p,2),size(p,2),0);
%         [KK, FF, B, ud, MM] = applyBCNullSpacesa(femodel.K,femodel.A,femodel.F,femodel.Q,femodel.G,femodel.H,femodel.R,femodel.M); % 0.7s
%         femodel.K0=femodel.K;
%         femodel.K=KK;
%         femodel.F=FF;
%         femodel.M=MM;
%         femodel.B=B;
%         [femodel.Nu,femodel.Or]=pdenullorth(femodel.H);
%         femodel.ud=ud;
%     elseif strcmp(method,'update')
%         femodel(1).F=createGlobalKFs(p,t,'F',fmatrix);
%     end
% end
end