function [ output_args ] = fkt_createHydroTimeSeries( NameOfBatch, i, discharge)
%CREATEFEEDTIMESERIES Summary of this function goes here
%   Detailed explanation goes here
format long e

%keyboard
filename = strcat('discharge_',NameOfBatch,'_',num2str(i),'.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%% Experiment Batch Name: %s \n', NameOfBatch);
fprintf(fileID,'%% Run Number:            %i \n', i);
fprintf(fileID,'%% Col1:          Col2:\n');
fprintf(fileID,'%% Time(h)        Discharge (m3/s) \n');

%Start with 0 feed
fprintf(fileID,'%i %e\n', 0, discharge);
fprintf(fileID,'%i %e\n', 2000000, discharge);
fclose(fileID);

% save as qb_feed.txt
copyfile(filename,'UserDefinedHydrograph.txt')
end

