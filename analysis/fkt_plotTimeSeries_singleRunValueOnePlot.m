function [ output_args ] = fkt_plotTimeSeries_singleRunValue(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
    TimeArrayHour, DataArray, PlotTitle, PlotYLabel, PlotXLabel,PlotNameStart,NameOfBatch,offlinedata_dir,sizeOfUsedArray)
%FKT_PLOTTIMESERIES1 Summary of this function goes here
%   Detailed explanation goes here

%%

%fig_visible = 'on';

fig = figure('visible',fig_visible);
set(fig, 'Position', [plot_x plot_y plot_width plot_height]);
if ~isnan(DataArray)
    c = jet(size(DataArray,1));
    for i=1:size(DataArray,1)
        x = TimeArrayHour(1:sizeOfUsedArray(i));
        y = DataArray(i,1:sizeOfUsedArray(i));
        plot(x,y,'Color',c(i,:));
        %plot(TimeArrayHour(1:sizeOfUsedArray(i)),DataArray(i,1:sizeOfUsedArray(i)),'.',c)
        hold on;
        %legendA = 'End of Run';
        
    end
end
title(PlotTitle);
%legend(legendA,'Location','northwest')
xlabel(PlotXLabel);
ylabel(PlotYLabel);
hold off;
% save figure
plotname = strcat(PlotNameStart,PlotTitle,'.png');
print(fig,plotname,'-dpng');
end

