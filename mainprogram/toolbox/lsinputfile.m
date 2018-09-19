function meshname=lsinputfile(a)
%datapath=['..' filesep 'datasave'];
filelist=dir(['toolbox' filesep 'modelandrheology' filesep '*.xlsx']);
for i=1:length(filelist)
    fprintf('%02.0f %20s\n',i,filelist(i).name)
end
prompt = ['input file for ' a];
n=input(prompt);
meshname=filelist(n).name;
% datafilelist=dir([datafoldername(n).folder filesep datafoldername(n).name filesep '*.mat']);
% Construct a questdlg with three options
% Include the desired Default answer
% load([dl(end).folder filesep dl(end).name])
end