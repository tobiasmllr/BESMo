function [ output_args ] = fkt_createFeedTimeSeries( RunConfig, i, RunParam)
%CREATEFEEDTIMESERIES Summary of this function goes here
%   Detailed explanation goes here
format long e

% Check for errors
assert(RunConfig.P_duration <= RunConfig.P_distanceArrayInMin(i), 'Error: RunConfig.P_duration must be shorter than RunConfig.P_distanceArrayInMin(i)')

TBefore_min	= (RunConfig.T_before * 60); % Initial Adjustment Time in minutes
TAfter_min	= (RunConfig.T_after * 60); % Adjustment Time in minutes after pulses stop
TSim_min    = (RunConfig.T_sim * 60); % Sim Time in minutes during which pulses occur

stepN           = RunParam.NPrint*RunParam.NtoPrint;              % total time steps number
param_dt        = RunParam.dt;
pulseNode_width = RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode);
pulseNode_dx    = RunParam.dx_array(RunParam.pulsedfeednode);

%% NEW BACKGROUND FEED RATE AS Fraction OF TOTAL FEED!
if RunConfig.MariaClaudiaPulses
    fprintf('Recreating feed of Maria and Claudias experiments!');
    
    %% Set up Arrays:
    time = (0 : stepN - 1) * param_dt;
    FeedArrayPulsed = zeros(1,length(time));
    FeedArrayConst  = zeros(1,length(time));
    
    %% constant feed rate of 300kg over 40h converted to MpS:
    pulseMag_kgPm = 300 / RunParam.RhoS / pulseNode_width / pulseNode_dx;
    backgroundFeedMpS = pulseMag_kgPm / (40 * 60 * 60);
    
    % precondition with constant feed for 6000 hours:
    expStartTime_sec = (TSim_min - (280 * 60))*60;
    
    startFeed=find(time<=0,1,'last');
    endFeed=find(time>=expStartTime_sec - param_dt,1,'first');
    
    FeedArrayConst(startFeed:endFeed) = backgroundFeedMpS;
    
    
    %% extract sequencing and save it to feed arrays
    phaseStart = 0;
    phaseLength = 40;
    
    sequencesN = size(RunConfig.MariaClaudiaSequence,1);
    for seq=1:sequencesN
        phaseEnd = phaseStart + phaseLength;
        switch(char(RunConfig.MariaClaudiaSequence(seq)))
            case '0'
                % no feed over 40 hours:
                startFeed   = find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed     = find(time>expStartTime_sec - param_dt + (phaseEnd * 60 * 60) - param_dt,1,'first');
                FeedArrayConst(startFeed:endFeed) = 0;
            case 'C'
                % constant feed over 40 hours:
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseEnd * 60 * 60) - param_dt,1,'first');
                FeedArrayConst(startFeed:endFeed) = backgroundFeedMpS;
            case '1'
                % one pulse:
                pulseMag_kgPm = 300 / RunParam.RhoS / pulseNode_width / pulseNode_dx;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
            case '2'
                % two pulses
                pulseMag_kgPm = 150 / RunParam.RhoS / pulseNode_width / pulseNode_dx;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
                
                phaseStart = phaseStart + 20;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
            case '4'
                % four pulses
                pulseMag_kgPm = 75 / RunParam.RhoS / pulseNode_width / pulseNode_dx;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
                
                phaseStart = phaseStart + 10;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
                
                phaseStart = phaseStart + 10;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
                
                phaseStart = phaseStart + 10;
                startFeed=find(time<=expStartTime_sec + (phaseStart * 60 * 60),1,'last');
                endFeed=find(time>expStartTime_sec - param_dt + (phaseStart * 60 * 60) + (RunConfig.P_duration * 60) - param_dt,1,'first');
                FeedArrayPulsed(startFeed:endFeed) = pulseMag_kgPm / (RunConfig.P_duration * 60);
            otherwise
                error('Dont know this phase type specified in RunConfig.MariaClaudiaSequence')
        end
        phaseStart = phaseEnd;
    end
else
    pulseMag_kgPm      = RunConfig.P_magnitudeArr(i) / RunParam.RhoS / pulseNode_width;
    totalFeed       = pulseMag_kgPm * RunConfig.P_Number(i);
    
    backgroundFeedMpS = RunParam.constantFeedFract * totalFeed / (TSim_min * 60);
    pulseMag_m_new = (1 - RunParam.constantFeedFract) * pulseMag_kgPm; 
    fprintf('Background feed: %.2f%% -> %.3f kg/hour \n',RunParam.constantFeedFract*100, backgroundFeedMpS * RunParam.RhoS * pulseNode_width * pulseNode_dx * 3600)
    fprintf('  Pulsed supply: %i events with %.2fkg each every %.2f hours\n',RunConfig.P_Number(i), pulseMag_m_new * RunParam.RhoS * pulseNode_width * pulseNode_dx, RunConfig.P_distanceArrayInMin(i) / 60)
    
    disp(strcat(datestr(now),' Starting to generate feed time series in fkt_createFeedTimeSeriesN'))
    %% Set up Arrays:
    time = uint32((0 : stepN - 1) * param_dt);
    FeedArrayPulsed = (double(zeros(1,length(time))));
    FeedArrayConst  = (double(zeros(1,length(time))));
    
    %% TIMING OF PULSES:
    TArrayStartSec  = uint32(round([(TBefore_min* 60) + 1:...
        (RunConfig.P_distanceArrayInMin(i)* 60):...
        (TBefore_min+(RunConfig.P_Number(i) - 1)*RunConfig.P_distanceArrayInMin(i)* 60) + 1])) ;     % Pulse starting times in seconds
    TArrayEndSec    = uint32(TArrayStartSec + (RunConfig.P_duration .* 60) - 1);  % Pulse ending times in seconds
    
    TArrayStartStep  = uint32(round([(TBefore_min* 60 / param_dt) + 1:...
        (RunConfig.P_distanceArrayInMin(i)* 60 / param_dt):...
        (TBefore_min+(RunConfig.P_Number(i) - 1)*RunConfig.P_distanceArrayInMin(i)* 60 / param_dt) + 1])) ;     % Pulse starting times in seconds
    TArrayEndStep    = uint32(TArrayStartStep + (RunConfig.P_duration .* 60 / param_dt) - 1);  % Pulse ending times in seconds
    
    if RunConfig.P_duration == RunConfig.P_distanceArrayInMin(i)
        % Constant Feed run!
        startFeed=find(time<=TArrayStartSec(1),1,'last');
        endFeed=find(time>=TArrayEndSec(end),1,'first');
        if isempty(endFeed)
            endFeed = length(time);
        end
        FeedArrayPulsed(startFeed:endFeed) = double((1 - RunParam.constantFeedFract) * totalFeed / (TSim_min * 60));
    else
        lengthFeedT = RunConfig.P_duration .* 60 / param_dt;
        pulseFeedRate = pulseMag_m_new / (RunConfig.P_duration .* 60);
        for p=1:(RunConfig.P_Number(i))
%                         % old slow:
%                         tic
%                         startFeed=find(time<=TArrayStartSec(p),1,'last');
%                         endFeed=find(time>=TArrayEndSec(p),1,'first');
%             
%                          lengthFeed1 = time(endFeed)-time(startFeed); % in seconds
%                         FeedArrayPulsed(startFeed:endFeed) = double(pulseMag_m_new / lengthFeed1);
%                         toc
%             
            %% new fast:
           FeedArrayPulsed(TArrayStartStep(p):TArrayEndStep(p))= deal(pulseFeedRate);
        end
    end
    
    assert(length(FeedArrayPulsed)==length(time),'problem with pulse timing!');
    
    % add constant feed fraction:
    constFeedPart = backgroundFeedMpS;
    
    %% CONSTANT FEED TO NODES:
    StartConst = find(time<=(TBefore_min * 60)); % timestep to start constant feed
    EndConst   = find(time>=(TBefore_min + TSim_min + TAfter_min) * 60 - param_dt); % timestep to stop constant feed
    
    FeedArrayConst(StartConst:EndConst) = (double(ones(1,EndConst-StartConst+1) .* constFeedPart)); % add constant feed rate to pulsed feed array
end
disp(strcat(datestr(now),' Finished generating feed time series in fkt_createFeedTimeSeries'))

%% SAVE Variables in .MAT INSTEAD

% Check for errors
if RunConfig.MariaClaudiaPulses
    feedAskedfor = 3 * 300 / RunParam.RhoS / pulseNode_width / pulseNode_dx;
    feedSimulated= sum(FeedArrayPulsed) * param_dt;
else
    feedAskedfor    = pulseMag_m_new * RunConfig.P_Number(i) / pulseNode_width / pulseNode_dx;
    feedSimulated   = (double(sum(FeedArrayPulsed) * param_dt)+double(sum(FeedArrayConst) * param_dt));
end
assert(abs((feedAskedfor)-(feedSimulated))<0.000001, 'Error: Total mass fed over sim does not add up: abs((pulseMag_m * P_Numberi)-(sum(FeedArray) * param_dt))<0.0000001')
% save feed time series and aux parameters
Runnumber = i;
P_Numberi = RunConfig.P_Number(i);
P_magnitudeArri = RunConfig.P_magnitudeArr(i);
NameOfBatch = RunConfig.NameOfBatch;

filename = strcat('qb_feed_',RunConfig.NameOfBatch,'_',num2str(i),'.mat');

save(filename, 'FeedArrayConst','FeedArrayPulsed', 'NameOfBatch', 'Runnumber','P_Numberi','P_magnitudeArri');
copyfile(filename,'qb_feed.mat')

end

