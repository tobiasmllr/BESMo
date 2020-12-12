function [ output_args ] = fkt_createGSD(i,NameOfBatch,GSD_Sizes,GSD_Pfi_initialbed,GSD_Pfi_feed)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

% Check for errors
assert(length(GSD_Sizes) == length(GSD_Pfi_initialbed) , 'Error: size GSD_Pfi_initialbed ~= GSD_sizes');
assert(length(GSD_Sizes) == length(GSD_Pfi_feed) , 'Error: size GSD_Pfi_feed ~= GSD_sizes');

filename = strcat('GSD_InitialSurf_',NameOfBatch,num2str(i),'.txt');
createGSDfile(filename, i, NameOfBatch, GSD_Sizes, GSD_Pfi_initialbed)
copyfile(filename,'GSD_InitialSurf.txt')

filename = strcat('GSD_InitialSubSurf_',NameOfBatch,num2str(i),'.txt');
createGSDfile(filename, i, NameOfBatch, GSD_Sizes, GSD_Pfi_initialbed)
copyfile(filename,'GSD_InitialSubSurf.txt')

filename = strcat('GSD_pbfeed_',NameOfBatch,num2str(i),'.txt');
createGSDfile(filename, i, NameOfBatch, GSD_Sizes, GSD_Pfi_feed)
copyfile(filename,'GSD_pbfeed.txt')
end

function [] = createGSDfile(filename, i, NameOfBatch, GSD_Sizes, GSD_Pfi)
fileID = fopen(filename,'w');
fprintf(fileID,'%% Experiment Batch Name: %s \n', NameOfBatch);
fprintf(fileID,'%% Run Number:            %i \n', i);
fprintf(fileID,'%% Col1:   Col2:\n');
fprintf(fileID,'%% Di (mm) Percent finer (%%) in ASCENDING ORDER\n');

for j=1:length(GSD_Sizes)
    fprintf(fileID,'%.3f\t%.3f\n', GSD_Sizes(j), GSD_Pfi(j));
end

fclose(fileID);
end

