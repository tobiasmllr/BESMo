function [ output_args ] = fkt_printSubplotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
    subplotHandle, subplotName)
%FKT_PRINTSUBPLOTTOEPS Summary of this function goes here
%   Detailed explanation goes here

if ~strcmp(version('-release'), '2013b')
    %% Plot single run figures
    %fig_visible = 'on';
    
    % QB DIV
    figSingle = figure('visible',fig_visible);
    %set(figSingle,'DefaultAxesFontSize',plot_fontsize)
    %set(figSingle,'DefaultTextFontSize',plot_fontsize)
    %set(figSingle,'Units','Inches');
    %set(figSingle, ...
    %    'Position', plot_position,...
    %    'PaperPosition', plot_paperposition, ...
    %    'PaperUnits', 'Inches', ...
    %    'PaperSize', plot_papersize)
    plotSingle = copyobj(subplotHandle,figSingle);
    set(plotSingle, 'Position', get(0, 'DefaultAxesPosition'));
    
    % omit title:
    set(get(gca,'Title'),'Visible','off');
    
    subplotName   = strcat(subplotName,'.tex');
    cleanfigure;
    matlab2tikz(subplotName,'showInfo', false,'strict',false);
    
    %print(figSingle,subplotName,'-depsc','-painters');
end
end

