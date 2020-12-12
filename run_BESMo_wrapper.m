function [] = run_BESMo_wrapper(nameconfigfile_batch, i, scratchDir, overwriteResults)
% This file will be executed by each worker individually

% if input number was interpreted as string, convert:
if ~isnumeric(i)
    i = str2num(i);
end

load(nameconfigfile_batch);

% check if this is one of the old configurations
if exist('RunConfig')
    Fnames = fieldnames(RunConfig);
    for field=1:size(Fnames,1)
        RunParam.(Fnames{field}) = RunConfig.(Fnames{field});
    end
end


if ~strcmp(RunParam.RunType,'LosPadres')
    [RunParam] = fkt_generateBatch(RunParam);
    
    if RunParam.continuePreviousConfig
        previous_run = strcat(RunParam.PreviousRunName,'_',int2str(i));
    else
        previous_run = '';
    end
else
    previous_run = '';
end

outputStyle = RunParam.outputStyle;
nc_params   = RunParam.nc_params;
run_mode    = RunParam.run_mode;

if run_mode.debug && run_mode.debugGenFilename
    timestring = dec2hex(round(now*1000000));
    disp(strcat('DEBUG MODE: Saving with appendix: ',timestring))
    current_run = timestring;
else
    current_run = strcat(RunParam.NameOfBatch,'_',int2str(i));
end

%% ========================================
% RUN MODEL
% ========================================
disp(strcat(' Preparing run:  ',current_run));
hostname = char( getHostName( java.net.InetAddress.getLocalHost ));
disp(strcat('hostname: ',hostname));
disp(strcat('scratch directory for temp storage: ',scratchDir));

% define output filename
resultfile = strcat(current_run,'.mat');

if ~exist(resultfile) || overwriteResults
    if run_mode.pickup && exist(strcat(current_run,'_varDump.mat'),'file')
        % Trying to pick up this configuration from an aborted state.
        % We will reload the old *_varDump.mat with all the feed/gsd/hydro
        % data.
        msg = strcat(datestr(now),' Found varDump.mat, picking up old state in run:  ',current_run);
        disp(msg);
    else
        if strcmp(RunParam.RunType, 'LosPadres')
            % CREATE GSD
            %fkt_createGSD_LosPadres(i,RunParam,RunParam);
            
            % CREATE FEED TIME SERIES
            %RunParam = fkt_createFeedTimeSeries_LosPadres(RunParam);
            
            % CREATE DISCHARGE TIME SERIES --- MOVED TO GENCONF
            %RunParam = fkt_createHydroTimeSeries_LosPadres(RunParam);
        else
            %% Prepare Pulse Run
            % CREATE GSD
            fkt_createGSD(i,RunParam.NameOfBatch,RunParam.GSD_sizes,RunParam.GSD_Pfi_initialbed,RunParam.GSD_Pfi_feed);
            
            % CREATE FEED TIME SERIES
            fkt_createFeedTimeSeries( RunParam, i, RunParam);
            
            % CREATE DISCHARGE TIME SERIES
            fkt_createHydroTimeSeries( RunParam.NameOfBatch, i, RunParam.Discharge);
        end
    end
    
    % CALL THE MODEL:
    msg = strcat(datestr(now),' Starting Single-run:  ',current_run);
    disp(msg);
    %fprintf(fdev, msg);
    if strcmp(RunParam.RunType, 'LosPadres')
        BESMo_LosPadres( RunParam, current_run, outputStyle, nc_params, run_mode, scratchDir, i)
    else
        BESMo( RunParam, current_run, previous_run, outputStyle, nc_params, run_mode, scratchDir)
    end
else
    msg = strcat(datestr(now),' Skipping run,  ',resultfile,' exists!');
    disp(msg);
    %fprintf(fdev, msg);
end

return
% RESULT: is a .mat file in the run folder. I load these in an other
% script to access the variables saved in there.
end
