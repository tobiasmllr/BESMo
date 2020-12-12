function [ output_args ] = fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
    plotHandle, plotName)
%FKT_PRINTSUBPLOTTOEPS Summary of this function goes here
%   Detailed explanation goes here

%% Plot single run figures
%fig_visible = 'on';

% QB DIV
%figSingle = figure('visible',fig_visible);
%set(figSingle,'DefaultAxesFontSize',plot_fontsize)
%set(figSingle,'DefaultTextFontSize',plot_fontsize)
%set(figSingle,'Units','Inches');
%set(figSingle, ...
%    'Position', plot_position,...
%    'PaperPosition', plot_paperposition, ...
%    'PaperUnits', 'Inches', ...
%    'PaperSize', plot_papersize)

%set(plotHandle, 'Position', get(0, 'DefaultAxesPosition'));

plotName   = strcat(plotName,'.eps');

    
% ommit title:
set(get(gca,'Title'),'Visible','off');
%set(gca,'Position',plot_position);
try
    cleanfigure;
catch
    warning('Add module matlab2tikz to path to allow cleanfigure command')
end
try
    print(plotHandle,plotName,'-depsc','-painters');
catch
    error('output failed, you might have to add export_fig to path https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig')
end

