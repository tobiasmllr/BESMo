function[DSize,DPercentFiner]=importGSD(filename,dataRow)
%Function to import the Surface GSD and convert it to a format which
%is used in the following calculations.
%File: column 1: Grain size in mm, column 2: Percent Finer in %
%Tobias July 11, 2013
format long e

%% load sediment sample
sedSampleInput = load(filename);

L = length(sedSampleInput);

%turn arround arrays:
DSize = sedSampleInput(1:L,1); % in mm
DPercentFiner = sedSampleInput(1:L,dataRow)/100;
end