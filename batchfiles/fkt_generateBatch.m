function [ RunParam ] = fkt_generateBatch( RunParam )
% GENERATE BATCH
% Uses RunConfig values to prepare batch runs
% RunConfig are parameters only used in the model set up,
% RunParam are used in the model itself
format long e

if ~RunParam.continuePreviousConfig
    if strcmp(RunParam.RunType,'A')
        RunParam = fkt_generateBatch_TypeA(RunParam);
    elseif strcmp(RunParam.RunType,'An')
        RunParam = fkt_generateBatch_TypeAn(RunParam);
    else
        RunParam = fkt_generateBatch_TypeB(RunParam);
        %else
        %    ME = MException('VerifyInput:OutOfBounds', ...
        %        strcat('Input Type not recognized: ',RunType));
        %    throw(ME);
    end
else
    RunParam.P_magnitudeArr       = RunParam.previous_P_magnitudeArr;
    RunParam.P_distanceArrayInMin = RunParam.previous_P_distanceArrayInMin;
    RunParam.P_Number             = round(RunParam.P_magnitude ./ RunParam.previous_P_magnitudeArr);
end
%% Store GSD
[ ~, ~, RunParam.GSD_sizes, RunParam.GSD_Pfi_feed ] = fkt_generateGSD(RunParam.GS_Classes,RunParam.GS_D50,RunParam.GS_Sigma, RunParam.GSDistributionType);
RunParam.GSD_Pfi_initialbed = RunParam.GSD_Pfi_feed;

end

