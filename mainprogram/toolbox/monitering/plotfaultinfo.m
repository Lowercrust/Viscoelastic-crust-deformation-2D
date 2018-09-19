function [faultstress,faultnodesdepth,varargout]=plotfaultinfo(model,stressyx,yq1,stresslimitplot,earthquakecount,coseismiccount,elapsedtime,varargin)
global modc
rho=modc.ther{1}.rho;
faultnodes=find(model.Mesh.Nodes(1,:)==0); % find nodes on the fault
faultnodesdepth=model.Mesh.Nodes(2,faultnodes);
faultstress=stressyx(faultnodes);
stresslowerlimitplot=stresslimitplot-5e6;

[faultdepth,plotdata]=findfaultdepth(model.Mesh,stressyx); % outputing plot data for stress
if length(varargin)==1
    u=varargin{1};
    faultslip=u(faultnodes);
    varargout{1}=faultslip;
    %faultstressdrop=stressdropyx(faultnodes);
    [faultnodesdepth,faultstress,faultslip]=sortdata(faultnodesdepth,faultstress,faultslip);
    plot(faultnodesdepth,faultstress,'r',yq1,stresslimitplot,'c',yq1,stresslowerlimitplot,'k',faultnodesdepth,-faultslip*1e8,'g',[faultdepth,faultdepth],[0,1e10],':');
    title(['EarthquakeNo.' num2str(earthquakecount) '@Time' num2str(elapsedtime/(365*3600*24)) '/EventNo.' num2str(coseismiccount)])
else
    plot(plotdata(:,1),plotdata(:,2),'r',yq1,stresslimitplot,'c',yq1,stresslowerlimitplot,'k')
    [faultnodesdepth,faultstress]=sortdata(faultnodesdepth,faultstress);
    title(['EarthquakeNo.' num2str(earthquakecount) 'Time' num2str(elapsedtime/(365*3600*24))])
end
if isempty(faultdepth)
    ylim([0,1e7]);
    xlim([0,1e7/(9.8*rho)]);
else
    ylim([0,rho*9.8*faultdepth*2*modc.mu+5e6]);
    xlim([0,faultdepth*2]);
end
drawnow
end