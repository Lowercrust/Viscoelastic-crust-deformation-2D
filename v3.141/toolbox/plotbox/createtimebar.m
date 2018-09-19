% input sorted timetable
function createtimebar(timetable,range,savefolder)
timetable=timetable(range(1):range(2))/(365*24*3600*1000000);
maxtime=timetable(end);
mkdir([savefolder filesep 'time'])
fig=figure('Position',[0 100 800 80]);
for i=1:length(timetable)
    barh(0.5,timetable(i));
    xlim([0 maxtime])
    ylim([0.2 0.7])
    set(gca,'ytick',[])
    xlabel('Time(Myr)')
%     axis tight
    saveas(fig,[savefolder filesep 'time' filesep 'time_' num2str(i,'%04.0f') '.png'])
end
end