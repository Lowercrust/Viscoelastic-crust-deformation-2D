function createmodelinfo(time)
global modc
modelinfo=dir('modelinfo.txt');
if isempty(modelinfo)
    fileID = fopen('modelinfo.txt','w');
    fprintf(fileID,'%15s %12s %12s %12s %12s %12s %12s %60s\r\n','Filename','Mesh nodes','Model width', 'Model depth', 'Velocity','Rheology','Mechanisms','Comments(60c)');
    % Rheology: WA: Wet anorthite , DA, Wet anorthite, WO:Wet olivine,
    % DO:Dry olivine....
else
    fileID = fopen(['modelinfo.txt'],'a');
end
prompt = 'Any comments or notes?(y or n)';
answer=input(prompt,'s');
if strcmp(answer,'y')
    prompt = 'please input the contents\n';
    contents = input(prompt,'s');
elseif strcmp(answer,'n')
    contents = '';
end
rhe='';
for i=1:length(modc.rhe)
    rhe=[rhe,modc.rhe{i}.Rheology];
end
fprintf(fileID,'%15s %12d %12d %12d %12.1f %12s %12s %60s\r\n',time,length(modc.Nodes),modc.mod.X,modc.mod.Y,modc.mod.v0,rhe,'TBA',contents);
fclose(fileID);