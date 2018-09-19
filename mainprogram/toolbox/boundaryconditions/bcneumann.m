% Neumann boundary condition for stress drop
function bcMatrix = bcneumann(problem,region,state)
global modc
%     scatter(region.x,region.y)
%     hold on
bcMatrix=-stressdropfun(region.y)/(modc.rhe{1}.G);
end

function F = stressdropfun(x)
% Stress drops to the stress lower limit
global stressyxi modc
stress=gather(interp1(stressyxi(:,1),stressyxi(:,2),x));
stressyxi=[];
if length(modc.ther)==1
    F = stress-(modc.ther{1}.rho*modc.g*x+modc.atm)*modc.mu;
elseif length(modc.ther)==2
    %not allow earthquake in mantle
    stresslimitmoho=(modc.ther{1}.rho*modc.g*modc.mod.CT+modc.atm)*modc.mu;
    if exist('modc.allowmantleearthquake','var')
        if length(x)==1
            if modc.allowmantleearthquake
                if x<modc.mod.CT
                    F = stress-(modc.ther{1}.rho*modc.g*(x)+modc.atm)*modc.mu;
                else
                    F= stress-((modc.ther{2}.rho*modc.g*(x-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho);
                end
            else
                if x<modc.mod.CT
                    F = stress-(modc.ther{1}.rho*modc.g*(x)+modc.atm)*modc.mu;
                else
                    F= 0;
                end
            end
        else
            %         disp(x);
            if x<modc.mod.CT
                F = stress-(modc.ther{1}.rho*modc.g.*(x)+modc.atm)*modc.mu;
            else
                F= stress-((modc.ther{2}.rho*modc.g.*(x-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho);
            end
        end
    else
        F=zeros(1,length(x));
        if exist('modc.allowmantleearthquake','var')
            if modc.allowmantleearthquake
               
                F(x<modc.mod.CT) = stress(x<modc.mod.CT) -(modc.ther{1}.rho*modc.g.*(x(x<modc.mod.CT))+modc.atm)*modc.mu;
                F(x>=modc.mod.CT)= stress(x>=modc.mod.CT)-((modc.ther{2}.rho*modc.g.*(x(x>=modc.mod.CT)-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho);
            else
  
                F(x<modc.mod.CT) = stress(x<modc.mod.CT) -(modc.ther{1}.rho*modc.g.*(x(x<modc.mod.CT))+modc.atm)*modc.mu;
                F(x>=modc.mod.CT)= 0;
            end
        else
            
            F(x<modc.mod.CT) = stress(x<modc.mod.CT) -(modc.ther{1}.rho*modc.g.*(x(x<modc.mod.CT))+modc.atm)*modc.mu;
            F(x>=modc.mod.CT)= stress(x>=modc.mod.CT)-((modc.ther{2}.rho*modc.g.*(x(x>=modc.mod.CT)-modc.mod.CT)+modc.atm)*modc.mu+stresslimitmoho);
        end
    end
    F(F<0)=0;
end
% F = interpolateSolution(resultstressyx,0,x)-rho*9.8*x*0.6;
% if x<faultdepth-100 && F>1e6
%     F=1e6;
% end
% if F<0
%     F=0;
% end
end