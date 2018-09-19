clear
close all
load('G:\result10241')
%% npointrect, stress boundary conditions on slip region, zero displacement below slip region
%     faultstress = interpolateSolution(createPDEResults(model,stressyx),xq1,yq1);
stresslimit=(rho*9.8*yq1*0.6+5e6)';
stresslowerlimit=(rho*9.8*yq1*0.6)';
faultstress = interpolateSolution(createPDEResults(model,stressyx),xq1,yq1); %wrong
%         stresslimit=(rho*9.8*yq1*0.6+5e6)';
%         stresslowerlimit=(rho*9.8*yq1*0.6)';
for jj=1:1
    gn=@npointrect;
    overlimit=(faultstress-stresslowerlimit);
    %         slip=find(overlimit>0);
    % [overlimit,sliplength]=cutNegativeToZero(overlimit);
    %         stressdroponfault=griddedInterpolant(yq1,overlimit);
    %         depth0=(slip(end)-1)*50;
    %     options = optimset('Display','iter');
    %         turn = fzero(@turningpoint,depth0);
    %     x1=0:3/(length(slip)-1):3;
    %     y1=1-erf(x1);
    overlimit(overlimit<0)=0;
    overlimit(isnan(overlimit))=0;
    [pointlist,nbs]=findsection(overlimit);
    %     overlimit(overlimit>0)=overlimit(overlimit>0).*y1';
    stressdroponfault=griddedInterpolant(yq1,overlimit);
    %slip region
    % pointlist=[0;turn;mod.Y];
    % nbs=length(pointlist)+2;
    % create model
    modeln = createpde;
    geometryFromEdges(modeln,gn);
    % apply boundary condition
    if rem(nbs,2)==1
        applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
        if nbs-4>2
            applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'u',0);
        end
        applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'g',@bcneumann,'Vectorized','on');
    else
        applyBoundaryCondition(modeln,'edge',[nbs-3,nbs-1],'u',0);
        if nbs-4>1
            applyBoundaryCondition(modeln,'edge',1:2:nbs-4,'u',0);
        end
        applyBoundaryCondition(modeln,'edge',2:2:nbs-4,'g',@bcneumann,'Vectorized','on');
    end
    
    % geometryFromEdges(modeln,g);
    %     pdegplot(modeln, 'edgeLabels', 'on');
    %         applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
    %         applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','on');
    % applyBoundaryCondition(modeln,'edge',3,'u',0);
    % applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','on');
    specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
    % generateMesh(modeln,'Hmax',100,'GeometricOrder','quadra','Jiggle','on','MesherVersion','R2013a');
    % resultn = solvepde(modeln);
    [u,mesh]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',50000,'ngen',inf); % coseismic slip
    %% change stressyx and stressyz to adapted Mesh
    stressyxc=meshchange(modeln,model,stressyx);
    stressyzc=meshchange(modeln,model,stressyz);
    %% result
    resultn=createPDEResults(modeln,u);
    stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake
    stressdropyz=rhe.G*resultn.YGradients;
    stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2);
    %         stressdropyxc=meshchange(model,modeln,stressdropyx); % change modeln's mesh to model's mesh
    %         stressdropyzc=meshchange(model,modeln,stressdropyz);
    %         stressdropc=meshchange(model,modeln,stressdrop);
    stressyxc=stressyxc-stressdropyx;
    stressyzc=stressyzc-stressdropyz;
    faultstress = interpolateSolution(createPDEResults(modeln,stressyxc),xq1,yq1); % shear stress tau_{yx}
    %         resultstressdroponfault = interpolateSolution(createPDEResults(model,stressdropc),xq1,yq1);
    %         resultstressdroponfaultyx = interpolateSolution(createPDEResults(modeln,stressdropyx),xq,yq);
    %         resultstressdroponfaultyz = interpolateSolution(createPDEResults(model,stressdropyzc),xq,yq);
    plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k')
    if max(faultstress)>0
        ylim([0,max(faultstress)*4]);
        xlim([0,max(faultstress)*4/(rho*9.8)]);
    end
    if isempty(overupperlimit(overupperlimit>0))
        break
    end
end
%         title(['earthquake' num2str(j) ',' num2str(jj)])
%         drawnow
%         anime2(j) = getframe(gcf);
%         maxfaultstress=max(faultstress);
%         minfaultstress=min(faultstress);
%         %         dstress=stress-stress0;
%         disp([maxfaultstress minfaultstress minetaeff dt/(365*24*3600) elapsedtime/(365*24*3600)])
%         overupperlimit=faultstress-stresslimit;

%         if isempty(overupperlimit(overupperlimit>0))
%             break
%         end

%         gn=@npointrect;
%         overlimit=(faultstress-stresslowerlimit);
%         slip=find(overlimit>0);
%
%         % [overlimit,sliplength]=cutNegativeToZero(overlimit);
%         %         stressdroponfault=griddedInterpolant(yq1,overlimit);
%         %         depth0=(slip(end)-1)*50;
%         %     options = optimset('Display','iter');
%         %         turn = fzero(@turningpoint,depth0);
%         %     x1=0:3/(length(slip)-1):3;
%         %     y1=1-erf(x1);
%         overlimit(overlimit<0)=0;
%         [pointlist,nbs]=findsection(overlimit);
%         %     overlimit(overlimit>0)=overlimit(overlimit>0).*y1';
%         stressdroponfault=griddedInterpolant(yq1,overlimit);
%         % slip region
%
%         %         pointlist=[0;turn;mod.Y];
%         % nbs=length(pointlist)+2;
%         % create model
%         modeln = createpde;
%         geometryFromEdges(modeln,gn);
%         % geometryFromEdges(modeln,g);
%         %     pdegplot(modeln, 'edgeLabels', 'on');
%         applyBoundaryCondition(modeln,'edge',[2,4],'u',0);
%         applyBoundaryCondition(modeln,'edge',1,'g',@bcneumann,'Vectorized','on');
%         % applyBoundaryCondition(modeln,'edge',2,'u',0);
%         % applyBoundaryCondition(modeln,'edge',4,'g',@bcneumann);
%         specifyCoefficients(modeln,'m',0,'d',0,'c',1,'a',0,'f',0);
%         [u,mesh]=generateadaptMesh(modeln,'MesherVersion','R2013a','maxt',10000,'ngen',inf); % coseismic slip
%         %% change stressyx and stressyz to adapted Mesh
%         stressyxc=meshchange(modeln,model,stressyx);
%         stressyzc=meshchange(modeln,model,stressyz);
%         %% result
%         resultn=createPDEResults(modeln,u);
%         stressdropyx=rhe.G*resultn.XGradients; % stressdrop after the earthquake
%         stressdropyz=rhe.G*resultn.YGradients;
%         stressdrop=sqrt(stressdropyx.^2+stressdropyz.^2);
% %         stressdropyxc=meshchange(model,modeln,stressdropyx); % change modeln's mesh to model's mesh
% %         stressdropyzc=meshchange(model,modeln,stressdropyz);
% %         stressdropc=meshchange(model,modeln,stressdrop);
%         stressyxc=stressyxc-stressdropyx;
%         stressyzc=stressyzc-stressdropyz;
%         faultstress = interpolateSolution(createPDEResults(modeln,stressyxc),xq1,yq1); % shear stress tau_{yx}
% %         resultstressdroponfault = interpolateSolution(createPDEResults(model,stressdropc),xq,yq);
% %         resultstressdroponfaultyx = interpolateSolution(createPDEResults(model,stressdropyxc),xq,yq);
% %         resultstressdroponfaultyz = interpolateSolution(createPDEResults(model,stressdropyzc),xq,yq);
%         plot(yq1,faultstress,'r',yq1,stresslimit,'c',yq1,stresslowerlimit,'k')
%         if max(faultstress)>0
%             ylim([0,max(faultstress)*5]);
%             xlim([0,max(faultstress)*5/(rho*9.8)]);
%         end
%         title(['earthquake' num2str(j) ',' num2str(jj)])
%         drawnow
%         anime2(j) = getframe(gcf);
%         maxfaultstress=max(faultstress);
%         minfaultstress=min(faultstress);
% %         dstress=stress-stress0;
%         disp([maxfaultstress minfaultstress minetaeff dt/(365*24*3600) elapsedtime/(365*24*3600)])
%         overupperlimit=faultstress-stresslimit;