function [ output_args ] = fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
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

set(plotHandle, 'Position', get(0, 'DefaultAxesPosition'));

% ommit title:
set(get(gca,'Title'),'Visible','off');

%set(gca,'Units','Inch');
%set(gca,'Position',plot_position);

plotName   = strcat(plotName,'.tex');
%     cleanfigure;
try
    matlab2tikz(plotName,'showInfo', false,'strict',false,'width','10cm','height','6cm' );
catch
    warning('Add module matlab2tikz to path to allow matlab tikz output for latex')
end
%print(figSingle,subplotName,'-depsc','-painters');
end

