function [ RunConfig ] = fkt_generateBatch_TypeA(RunConfig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

%% Calculate parameters
% Pulse TIMES in experiment
P_distanceArrayInMinAll = [...
    RunConfig.P_distanceStart:...
    RunConfig.P_distanceIncrease:...
    (RunConfig.T_sim * 60)];

%Number of pulses:
% Only keep integer numbered NPulse
j = 1;
for i = 1:length(P_distanceArrayInMinAll)
    P_Numberi = (RunConfig.T_sim * 60) ./ P_distanceArrayInMinAll(i);
    P_Numberi = P_Numberi(P_Numberi >= RunConfig.P_NMin);
    if ceil(P_Numberi) == floor(P_Numberi)
        P_Number(j) = P_Numberi;
        P_distanceArrayInMin(j) = P_distanceArrayInMinAll(i);
        j = j + 1;
    end
end

%Magnitude of each pulse, minimum P_NMin pulses
P_magnitudeArr = RunConfig.P_magnitude ./ P_Number;

% save parameters into structure:
RunConfig.P_distanceArrayInMin = P_distanceArrayInMin;
RunConfig.P_Number = P_Number;
RunConfig.P_magnitudeArr = P_magnitudeArr;

%save(NameOfConfig)
end

