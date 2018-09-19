% plot stress
function pdeplotsavedata(data,model,field,time,num,period,figurefolder,stresslimit)
switch field
    % if strcmp(field,'tyx')
    case 'tyx'
        h=pdeplot(model,'xydata',data,'colormap','parula');
%         colormap(parula)
        %     colorlimit=varargin{1};
        caxis([stresslimit(1,1) stresslimit(1,2)])
        titlestring='$\tau_{yx}$';
        % elseif strcmp(field,'tyz')
    case 'tyz'
        h=pdeplot(model,'xydata',data,'colormap','parula');
%         colormap(parula)
        %     colorlimit=varargin{1};
        caxis([stresslimit(2,1) stresslimit(2,2)])
        titlestring='$\tau_{yz}$';
    case 'result'
        if  strcmp(period,'Interseismic')
            titlestring='Interseiemic velocity';
            h=pdeplot(model,'xydata',(data)*365*24*3600*1000);
            field='result';
        elseif  ~strcmp(period,'Interseismic')
            titlestring='Coseismic displacement';
            h=pdeplot(model,'xydata',data);
            field='result';
        end
    case 'etaeff'
        %     Fdisl=varargin{1};
        %[etaeff,esyxv,esyzv]=stress2visco(Fdisl,data.tyx,data.tyz);
        %     data=removeinf(data);
        %     if strcmp(field,'etaeff')
        titlestring='$\eta_{eff}$';
        h=pdeplot(model,'xydata',removeinf(log10(abs(data))));
        colormap(flipud(jet))
        caxis([19 25])
    case 'esyxv'
        %     elseif strcmp(field,'esyxv')
        h=pdeplot(model,'xydata',removeinf(log10(data)));
        titlestring='$\dot{\epsilon}_{yx}$';
        caxis([-20 -14])
    case 'esyzv'
        %     elseif strcmp(field,'esyzv')
        h=pdeplot(model,'xydata',removeinf(log10(abs(data))));
        titlestring='$\dot{\epsilon}_{yz}$';
        caxis([-20 -14])
end
% h.Position=[100 100 600 800];
%     if ~isempty(varargin)
%         colorlimit=varargin{1};
%         if size(colorlimit,2)==2
%             caxis([colorlimit(1) colorlimit(2)])
%         end
%     end
axis ij %
axis equal tight
title([titlestring ' at t=' num2str(round(time))],'interpreter','latex') %

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% savefilename=[figurefolder filesep field filesep period(1:2) filesep num2str(num,'%06.0f') '.png'];
% print(savefilename,'-dpng','-r0')
% drawnow
saveas(gcf,[figurefolder filesep field filesep period(1:2) filesep num2str(num,'%06.0f') '.png']);
end
