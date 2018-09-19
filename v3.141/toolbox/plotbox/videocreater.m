function videocreater(figpath,videoname,framerate)
v = VideoWriter([figpath filesep videoname]);
v.FrameRate = framerate;
open(v);
fl=dir([figpath filesep '*.png']);
disp(['convert ' num2str(length(fl)) ' figures into video with frame rate of ' num2str(framerate)]);
progressbar;
for i=1:length(fl)
    progressbar(i/length(fl))
    frame=imread([fl(i).folder filesep fl(i).name]);
    writeVideo(v,frame);
end
close(v);
end