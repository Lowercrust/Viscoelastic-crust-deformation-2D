function meshname=lsmeshfile
%datapath=['..' filesep 'datasave'];
meshnamelist=dir(['toolbox' filesep 'mesh' filesep '*.mat']);
for i=1:length(meshnamelist)
    fprintf('%02.0f %20s\n',i,meshnamelist(i).name)
end
prompt = 'Which mesh?';
n=input(prompt);
meshname=meshnamelist(n).name;
% datafilelist=dir([datafoldername(n).folder filesep datafoldername(n).name filesep '*.mat']);
% Construct a questdlg with three options
% Include the desired Default answer
% load([dl(end).folder filesep dl(end).name])
end