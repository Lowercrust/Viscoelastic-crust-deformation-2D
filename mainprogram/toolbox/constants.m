classdef constants
    properties
        mod; % model geometry
        rhe; % model rheology
        ther;% model thermal properties
        g,R,v0,atm,mu,fs,fw;
        % gravitational acceleration, gas constant, boundary
        % veolocity, atmospheric pressure, coefficient of friction, fault
        % strength, fault width
        Nodes,ms; % nodes and minimum grid size
        allowmantleearthquake;
        considerplastic
        % minimum grid size is only exist when creating new mesh
    end
    methods
        function obj=constants()
            global mod
            [~,~,model]=xlsread(lsinputfile('Model Properties'));
            %             [~,~,thermal]=xlsread(lsinputfile('Thermal Properties'));
            
            mod=cell2table(model(:,2)');
            mod.Properties.VariableNames=model(:,1)';
            mod=table2struct(mod);
            obj.mod=mod;
            
%             ther=cell2table(thermal(:,2)');
%             ther.Properties.VariableNames=thermal(:,1)';
%             obj.ther=table2struct(ther);
            n=input('(1) single layer crust \n(2) crust-mantle two layer\n');
            mu=input('fault frictonal coefficent');
            pn=input('(0)elastic\n(1)elastic-plastic with strain hardening\n(2)elastic-perfectly plastic\n');
            if n==2
                prompt = 'allow earthquakes in mantle(y or n)';
                answer=input(prompt,'s');
                if strcmp(answer,'y')
                    obj.allowmantleearthquake=true;
                elseif strcmp(answer,'n')
                    obj.allowmantleearthquake=false;
                end
            end
            obj.rhe=cell(n,1);
            obj.ther=cell(n,1);
            for i=1:n
                [~,~,rheprop]=xlsread(lsinputfile(['File for rheology ' num2str(i)]));
                rhe0=cell2table(rheprop(:,2)');
                rhe0.Properties.VariableNames=rheprop(:,1)';
                obj.rhe{i}=table2struct(rhe0);
                [~,~,therprop]=xlsread(lsinputfile(['File for thermal properties ' num2str(i)]));
                ther0=cell2table(therprop(:,2)');
                ther0.Properties.VariableNames=therprop(:,1)';
                obj.ther{i}=table2struct(ther0);
            end
            obj.considerplastic=pn;
            obj.R=8.3144598;%gas constant [J K^-1 mol^-1]
%             obj.rho=2800;%density[kg m^-3]\
            obj.g=9.8;
            obj.mu=mu;
            obj.fs=5e6;
            obj.fw=0.5;
            obj.atm=101325; % [pa]
            obj.v0=mod.v0/(365*3600*24*1000);%[m/s] boundary velocity
            prompt = 'create new mesh?(y or n)';
            answer=input(prompt,'s');
            if strcmp(answer,'y')
                %                 generateMesh2d(obj,'linerefinement');
                prompt = 'input minimum mesh size';
                obj.ms=input(prompt);
                nodes=generateNodes2d(obj);
                obj.Nodes=nodes;
            elseif strcmp(answer,'n')
                load(lsmeshfile);
                if exist('msh','var')
                    obj.Nodes=msh.Nodes;
                elseif exist('nodes','var')
                    obj.Nodes=nodes;
                elseif exist('vert','var')
                    disp('post processing')
                    p1=vert';
                    p1=deleteoutside(p1,[0,obj.mod.X,0,obj.mod.Y]);
                    [p1,~]=addboundarynodes(p1,obj.fw);
                    obj.Nodes=p1;
                end
            end
        end
    end
end

function nodes=generateNodes2d(obj)
global xs ms 
xs=obj.fw;
ms=obj.ms;
optr.kind='delaunay';
optr.rho2=1.1;
% optr.siz1=1;
% optr.siz2=10;
optr.disp=1;
opts.disp=1;
X=obj.mod.X;
Y=obj.mod.Y;
CT=obj.mod.CT;

node=[xs,0;X,0;X,CT;xs,CT;xs,Y;X,Y];
% edge=[1,5;5,6;6,2;2,1;1,4;4,3;3,2;2,1];
edge=[4,1;1,2;2,3;3,4;5,4;4,3;3,6;6,5];
part=cell(1,2);
part{1}=[1,2,3,4];
part{2}=[5,6,7,8];
%%
time=now;
% hfun1=@testhfun1;
[node,edge,~,~] = refine2(node,edge,part,optr,200);
% node=modifye5nodes(node);
% refine mesh on horizonal boundary
hfun=@boundaryrefinex;
[vert,etri,~,~] = refine2(node,edge,[],optr,hfun);
% refine mesh on vertical boundary
hfun=@boundaryrefiney;
[vert,etri,tria,tnum] = refine2(vert,etri,[],optr,hfun);
% refine mesh on moho
if ~isempty(obj.mod.CT)
    hfun=@boundaryrefinemoho;
    [vert,etri,tria,tnum] = refine2(vert,etri,[],optr,hfun);
end
time=(now-time)*24*3600;
disp(['mesh refined\elapsed time' num2str(time)])

%%
% disp('smoothing mesh')
% time=now;
% [vert,~,~,~] = smooth2(vert,etri,tria,tnum,opts);
% time=(now-time)*24*3600;
% disp(['mesh smoothed\elapsed time' num2str(time)])

%%
disp('post processing')
p1=vert';
p1=deleteoutside(p1,[0,obj.mod.X,0,obj.mod.Y]);
[p1,~]=addboundarynodes(p1,xs);
nodes=p1;
boundarynodeindex=find(p1(1,:)==0);
%% save mesh
formatOut = 'yymmddHHMMSS';
time=datestr(now,formatOut); % msh generated time.
meshpath='./toolbox/mesh/';
save([meshpath time '.mat'],'nodes')
%% output mesh info
meshlist=dir([meshpath 'meshinfo.txt']);
if isempty(meshlist)
    fileID = fopen([meshpath 'meshinfo.txt'],'w');
    fprintf(fileID,'%15s %12s %12s %12s %12s \r\n','filename','nodenumber', 'boundarynode','modelwidth','modeldepth');
else
    fileID = fopen([meshpath 'meshinfo.txt'],'a');
end
fprintf(fileID,'%15s %12d %12d %12d %12d\r\n',time,length(p1),length(boundarynodeindex),obj.mod.X,obj.mod.Y);
fclose(fileID);
%fprintf(fileID,fmt,[55 55 55 55]);
end

function p1=deleteoutside(p1,range)
poutsideindexx=find(p1(1,:)<range(1) | p1(1,:)>range(2));
poutsideindexy=find(p1(2,:)<range(3) | p1(2,:)>range(4));
dpoint=unique([poutsideindexx,poutsideindexy]);
p1(:,dpoint)=[];
end

function hfun=boundaryrefinex(test)
global ms mod
hmin=ms;
hmax=500;
hfun=ones(size(test,1),1)*hmax;
hfun(test(:,2)<=100)=(test(test(:,2)<=100,2))*(hmax-hmin)/100+hmin;% surface
hfun(test(:,2)>=mod.Y-100)=(mod.Y-(test(test(:,2)>=mod.Y-100,2)))*(hmax-hmin)/100+hmin;% bottom
end

function hfun=boundaryrefiney(test)
global xs ms mod
hmin=ms;
hmax=500;
hfun=ones(size(test,1),1)*hmax;
hfun(test(:,1)<=100)=(test(test(:,1)<=100,1)-xs)*(hmax-hmin)/100+hmin; % fault
hfun(test(:,1)>=mod.X-100)=(mod.X-(test(test(:,1)>=mod.X-100,1)))*(hmax-hmin)/100+hmin;% farfield
end

function hfun=boundaryrefinemoho(test)
global ms mod
hmin=ms;
hmax=500;
hfun=ones(size(test,1),1)*hmax;
hfun(test(:,2)<mod.CT+50)=(-mod.CT+(test(test(:,2)<mod.CT+50,2)))*(hmax-hmin)/50+hmin;% moho+50
hfun(test(:,2)>=mod.CT-50)=(mod.CT-(test(test(:,2)>=mod.CT-50,2)))*(hmax-hmin)/50+hmin;% moho-50
end



% copy nodes on E5 to boundary
function [p1,newnodesnum]=addboundarynodes(p1,xs)
% global xs
E5nodes=find(single(p1(1,:))==single(xs));
boundarynodes=find(p1(1,:)==0,1);
if isempty(boundarynodes)
    newnodes=[zeros(1,length(E5nodes));p1(2,E5nodes)];
    p1=[p1,newnodes];
    newnodesnum=length(newnodes);
    disp([num2str(newnodesnum) ' nodes have been added']);
else
    disp('Boundary nodes exist, no need to add')
    newnodesnum=0;
end
end
