function dispcalculationinfo(calinfo)
if isempty(calinfo.dt) || calinfo.dt==0
    if calinfo.count==1
        fprintf(1,'%10s %10s %10s %10s\n','count','timeforcal','faultdepth','maxslip');
    end
    fprintf(+1, '%10.0f %10.2f %10.2f %10.4f \n',...
        [calinfo.count,calinfo.timeforcal,calinfo.faultdepth,calinfo.maxslip]);
else
    if isempty(calinfo.dT)
        if calinfo.count==1
            fprintf(1,'%10s %10s %10s %10s %10s\n','count','timeforcal','time','dt','istime');
        end
        fprintf(+1, '%10.0f %10.2f %10.2f %10.2f %10.2f\n',...
            [calinfo.count,calinfo.timeforcal,calinfo.elapsedtime/(365*3600*24),calinfo.dt/(365*3600*24),calinfo.intertime/(365*3600*24)]) ;
    else
        if calinfo.count==1
            fprintf(1,'%10s %10s %10s %10s %10s %10s %10s\n','count','timeforcal','time','dt','istime','max(dT)','max(dT0)');
        end
        fprintf(+1, '%10.0f %10.2f %10.2f %10.2f %10.2f %10.2e %10.2e \n',...
            [calinfo.count,calinfo.timeforcal,calinfo.elapsedtime/(365*3600*24),calinfo.dt/(365*3600*24),calinfo.intertime/(365*3600*24),calinfo.dT,calinfo.dT0]) ;
    end
end
end