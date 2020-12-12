function [ RunConfig ] = fkt_generateBatch_TypeB(RunConfig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

%% Calculate parameters
% Pulse TIMES in experiment
P_distanceArrayInMinAll = [...
    RunConfig.P_distanceStart:...
    RunConfig.P_distanceIncrease:...
    (RunConfig.T_sim * 60)];

% Keep all pulses:
P_distanceArrayInMin = P_distanceArrayInMinAll;
for i = 1:length(P_distanceArrayInMinAll)
    P_Numberi(i) = (RunConfig.T_sim * 60) ./ P_distanceArrayInMinAll(i);
end
P_Number = P_Numberi(P_Numberi >= RunConfig.P_NMin);

% In case of constant magnitude:
% Set magnitude array to the size of NPulse and save the same magnitude for
% each pulse
P_magnitudeArr = P_Number;
P_magnitudeArr(:) = RunConfig.P_magnitude;

% save parameters into structure:
RunConfig.P_distanceArrayInMin = P_distanceArrayInMin;
RunConfig.P_Number = P_Number;
RunConfig.P_magnitudeArr = P_magnitudeArr;

%save(NameOfConfig)
end

