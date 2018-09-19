function exitflag=isinterseismic(interseismiccount,faultslip,faultstress,faultnodesdepth,stresslimit)
% global rho faultstrength faultdepth
global modc faultdepth
% disp([min(faultstress), max(faultstress)])
faultindex=find(faultnodesdepth<faultdepth);
fend=faultindex(end);
faultstress(max(1,fend-4000):fend)=smooth(faultstress(max(1,fend-4000):fend),1001);
if exist('modc.allowmantleearthquake','var')
    if modc.allowmantleearthquake
        index=find(faultstress'-stresslimit>0);
    else
        index=find(faultstress(faultnodesdepth<modc.mod.CT)'-stresslimit(faultnodesdepth<modc.mod.CT)>0);
    end
else
    index=find(faultstress'-stresslimit>0);
end
% index=find(faultstress'-stresslimit>0);
if isempty(index) %|| min(faultstress)<1e4 % interseismic
    exitflag=true; % when faultstress is smaller than stresslimit
elseif max(abs(faultslip))<1e-4 && interseismiccount==0
    exitflag=true; % when faultslip is very small, coseismic period ends
else
    exitflag=false; % when fault stress is larger then stresslimit
end

if max(abs(faultslip))>=1e-4
    exitflag=false;
end
end