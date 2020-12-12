function [  ] = fkt_plotTimeSeries_singleRunValue(i, TimeArray, DataArray, ...
    plotXmin, plotXmax, plotYmin, plotYmax, ...
    PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray)
%FKT_PLOTTIMESERIES1 Summary of this function goes here
%   Detailed explanation goes here

%% run
plot(TimeArray(:),DataArray(:),'-');

if(plotYmin<plotYmax)
    axis([plotXmin plotXmax plotYmin plotYmax]);
else
    warning('plotYmin >= plotYmax');
end
hold on;
%legendA = 'End of Run';

%title(PlotTitle)
%legend(legendA,'Location','northwest')
xlabel(PlotXLabel);
ylabel(PlotYLabel);
end

