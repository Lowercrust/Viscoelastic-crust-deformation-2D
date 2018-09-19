gn=@npointrect;
% find fault depth
overlimit=faultstress-stresslowerlimit;
% if isempty(pointlist)
[pointlist,nbs]=findsection(overlimit);
% else
%     pointlist0=pointlist(2:end-1);
%     [pointlist,~]=findsection(overlimit);
%     pointlist=sort([pointlist0 pointlist]);
%     nbs=length(pointlist)+2;
%     disp(pointlist)
% end
% pointlist=[0 35000];
% nbs=4;
modeln = createpde;
geometryFromEdges(modeln,gn);
% pdegplot(modeln, 'edgeLabels', 'on');
% apply boundary condition
if nbs==4
    applyBoundaryCondition(modeln,'edge',3,'u',0);
    applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','off');
else
    if rem(nbs,2)==1
        applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
        if nbs-4>2
            applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'u',0);
        end
        applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'g',@bcneumann,'Vectorized','off');
    elseif rem(nbs,2)==0
        applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
        if nbs-4>1
            applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'u',0);
        end
        applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'g',@bcneumann,'Vectorized','off');
    end
end
specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
% [p1,e1,t1]=meshToPet(model.Mesh);
% calculate displacement and stressdrop
[u,mesh]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',50000,'ngen',inf); % coseismic slip
resultn=createPDEResults(modeln,u);
stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake yx
stressdropyz=rhe.G*resultn.YGradients; % yz
stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2); % maximum shear stress
% interpolate stress to adaptmesh
tic
stressyx=meshchange(modeln,model,stressyx);
toc
stressyz=meshchange(modeln,model,stressyz);
toc
model=modeln;
% update stressfield
stressyx=stressyx-stressdropyx; 
stressyz=stressyz-stressdropyz;
stress=sqrt(stressyx.^2+stressyz.^2);
% plot fault stress
xq1 = linspace(0,0,mod.Y+1);
yq1 = linspace(0,mod.Y,mod.Y+1);
resultstressyx=createPDEResults(modeln,stressyx);
faultstress = interpolateSolution(resultstressyx,xq1,yq1); 
stresslowerlimit=(rho*9.8*yq1*0.6)';
stresslimit=(rho*9.8*yq1*0.6+1e6)';
faultslip=interpolateSolution(resultn,xq1,yq1);
plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k',yq1,-faultslip*1e9,'g')
if max(faultstress)>0
    ylim([0,max(faultstress)*2]);
    xlim([0,max(faultstress)*2/(rho*9.8)]);
end