function [StructPulses] = write_pulseSize_struct(nameconfigfile_batch, i, StructPulses)
% Write a table showing the pulse sizes for the different simulations
load(nameconfigfile_batch);
%RunParam = fkt_generateBatchConstants();
[RunConfig] = fkt_generateBatch(RunConfig);

if RunConfig.continuePreviousConfig
    StructPulses(i).Run = i;
    StructPulses(i).PulseNumber = round(RunConfig.P_magnitude ./ RunConfig.previous_P_magnitudeArr(i));
    StructPulses(i).PulseMagnit = RunConfig.previous_P_magnitudeArr(i);
    StructPulses(i).PulsePeriod = RunConfig.previous_P_distanceArrayInMin(i) ./ 60;
else
    pulseMag_m_new = (1 - RunParam.constantFeedFract) * RunConfig.P_magnitudeArr(i) / RunParam.RhoS / RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode) / RunParam.dx_array(RunParam.pulsedfeednode);
    pulseMag    = pulseMag_m_new * RunParam.RhoS * RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode) * RunParam.dx_array(RunParam.pulsedfeednode);
    pulsePeriod = RunConfig.P_distanceArrayInMin(i) / 60;
    %fprintf('  Pulsed supply: %i events with %.2fkg each every %.2f hours\n',RunConfig.P_Number(i), pulseMag, pulsePeriod)
    
    % append to table:
    StructPulses(i).Run = i;
    StructPulses(i).PulseNumber = RunConfig.P_Number(i);
    StructPulses(i).PulseMagnit = pulseMag;
    StructPulses(i).PulsePeriod = pulsePeriod;
end
end