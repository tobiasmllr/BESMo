function [ ] = fkt_analyzeRuns_importVar( current_run, batchDir, RunConfig, RunParam, WorkspaceFile, goThroughSubsurfaceLayers,StartAnalysisHour)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

nGSD = length(unique(RunConfig.GSD_sizes));

TBefore_h = RunConfig.T_before;

TExp_h = RunConfig.T_sim;
TstepDt = RunParam.dt;
ToutputFreq = RunParam.NtoPrint;
TPulseDelta_minArr = RunConfig.P_distanceArrayInMin;
NPulse = RunConfig.P_Number;
%keyboard
%% PULSE RUN PARAMETERS:
TPulseDelta_hArr = TPulseDelta_minArr / 60;

%% BATCH RUN PARAMETERS
% option to start analysis later in the run:
assert(StartAnalysisHour<TExp_h);

feed_starting_store     = 1 + (TBefore_h + StartAnalysisHour) * 3600 / TstepDt / ToutputFreq;
feed_starting_store_real= 1 + TBefore_h * 3600 / TstepDt / ToutputFreq;
feed_ending_store     = feed_starting_store_real + (TExp_h * 3600 / TstepDt / ToutputFreq);

feed_starting_step   = 1 + (TBefore_h + StartAnalysisHour) * 3600 / TstepDt;
feed_ending_step     = feed_starting_step - 1 + ((TExp_h - StartAnalysisHour) * 3600 / TstepDt);

feed_starting_sec   = (TBefore_h + StartAnalysisHour) * 3600;
feed_ending_sec     = feed_starting_sec + ((TExp_h - StartAnalysisHour) * 3600);

% INDEX FOR VARIABLES:
%run_index_tru = run_starting:run_overstep:run_ending;
arrSize = [length(NPulse) feed_ending_store - feed_starting_store];

%PREALLOCATION

%PREALLOCATE VAR
volume = NaN(arrSize(1),arrSize(2));
volumeDiff = NaN(arrSize(1),arrSize(2)-1);
volumeDiffCumsum = NaN(arrSize(1),arrSize(2)-1);
volumeDiffEnd = NaN(arrSize(1),1);
volumeDiffMax = NaN(arrSize(1),1);
volumeDiffMin = NaN(arrSize(1),1);
volumeDiffMean = NaN(arrSize(1),1);
% SLOPE (omitting first and last section
%slopeCells = NaN(arrSize(1),arrSize(2));
slopeDiffEnd = NaN(arrSize(1),1);
slopeDiffMax = NaN(arrSize(1),1);
slopeDiffMin = NaN(arrSize(1),1);
slopeDiffMean = NaN(arrSize(1),1);

% GRAIN SIZE DISTIRBUTIONS
GSD_Sizes       = NaN(arrSize(1),nGSD);
GSD_initSurfFsi = NaN(arrSize(1),nGSD);          % GSD init surf:
GSD_feedFsi     = NaN(arrSize(1),nGSD);          % GSD Feed:
%GSD Time series:
% single values:
GSD_TS_90       = NaN(arrSize(1),arrSize(2)); % GSD D90
GSD_TS_Dsg      = NaN(arrSize(1),arrSize(2)); % GSD Dsg
GSD_TS_Dssg     = NaN(arrSize(1),arrSize(2)); % GSD Dssg
% distributions:
GSD_TS_SurfFsi  = NaN(arrSize(1),nGSD,arrSize(2));
GSD_TS_TransPbi = NaN(arrSize(1),nGSD-1,arrSize(2));

% FEED Volume
%feedVol_m3_pulse = NaN(arrSize(1),1);
%feedVol_m3_const = NaN(arrSize(1),1);
%feedVol_kg_pulse = NaN(arrSize(1),1);
%feedVol_kg_const = NaN(arrSize(1),1);

% Time when the next pulse would hit if the simulation time was longer:
TPulseDistance_store = NaN(arrSize(1),1);
TPulseNextWouldHit_store = NaN(arrSize(1),1);
sizeOfUsedArray = NaN(arrSize(1),1);
% Qb Sediment Transport at 1 nodes from flume exit
qb_out = NaN(arrSize(1),arrSize(2));
% qb_in_const_step = NaN(arrSize(1),feed_ending_step-feed_starting_step+1);
% qb_in_pulse_step = NaN(arrSize(1),feed_ending_step-feed_starting_step+1);

qb_in_const_store = NaN(arrSize(1),arrSize(2));
qb_in_pulse_store = NaN(arrSize(1),arrSize(2));

qb_mean = NaN(arrSize(1),1);
%qb_div = NaN(arrSize(1),arrSize(2));
qb_divcum = NaN(arrSize(1),arrSize(2));
qb_divdiff = NaN(arrSize(1),arrSize(2)-1);
qb_output_m = NaN(arrSize(1),arrSize(2));
%qb_output_m3 = NaN(arrSize(1),arrSize(2));
%qb_output_kg = NaN(arrSize(1),arrSize(2));
% Average Qb over pulse-period
qb_mean_pulse = NaN(arrSize(1),max(NPulse));
pulse_timing_sec = NaN(arrSize(1),max(NPulse));
%percentual QbDiffCum
qb_divcumPerc = NaN(arrSize(1),arrSize(2));
%slope
slope       = NaN(arrSize(1),arrSize(2));
slope_mean  = NaN(arrSize(1),1);
% ustar
ustar       = NaN(arrSize(1),arrSize(2));
ustar_mean  = NaN(arrSize(1),1);

% pulses
nPulseBegins_store  = NaN(arrSize(1),max(NPulse));
nPulseEnds_store    = zeros(arrSize(1),max(NPulse));
% loop
%% batch run loop

for i=1:length(NPulse)
    close all
    
    runDir = strcat(batchDir,'/',current_run{1,i});
    nc_filename     = strcat(batchDir,'/',current_run{1,i},'_store.nc');
    mat_filename    = strcat(batchDir,'/',current_run{1,i},'.mat');
    if (~exist(nc_filename,'file') || ~exist(mat_filename,'file'))
        nc_filename  = strcat(runDir,'/',current_run{1,i},'_store.nc');
        mat_filename = strcat(runDir,'/',current_run{1,i},'.mat');
        if (~exist(nc_filename,'file') || ~exist(mat_filename,'file'))
            msg = strcat('Error: File not found: ',nc_filename);
            warning(msg)
        end
    end
    fprintf(strcat('import of analysis NetCDF-file: ',nc_filename,'... \n'))
    
    %fprintf(strcat('Extracting Variables from run ',num2str(i),'... \n'))
    if (exist(nc_filename,'file') && exist(mat_filename,'file'))
        
        % load non-temporal arrays from mat-file (no historical data)
        % use regexp to exclude the variables that start with the listed
        % names, as these contain data for every timestep
        rundata = load(mat_filename, '-regexp', '^(?!Courant|hydro_TArray|subtimestepsArray|overThreshold|time)...');
        
        % import the historical data into rundata-structure by importing
        % NetCDF data
        [ rundata ] = netcdf_import_all( nc_filename, rundata, goThroughSubsurfaceLayers );
                
        % FOCUS NODE:
        disp('focusArray now taken as rundata.nodes_loc(RunParam.pulsedfeednode:end)')
        focusArray = RunParam.pulsedfeednode+1:length(rundata.nodes_loc)-2;
        %focusArray = 
        
        focusExportNode = length(rundata.nodes_loc) - 1; % output node
        % ELEVATION (in m)
        elevCells = rundata.store_etab;
        
        % VOLUME (in m3)
        %keyboard
        volume(i,:)  = sum(elevCells(:,feed_starting_store:feed_ending_store-1) .* repmat(rundata.param_reachbankfull .* rundata.dx_array,1,feed_ending_store-feed_starting_store));
        volumeDiff(i,:) = diff(volume(i,:));
        volumeDiffCumsum(i,:) = cumsum(volumeDiff(i,:));
        volumeDiffEnd(i) = (volume(i,end) / volume(i,1));
        volumeDiffMax(i) = (max(volume(i,:))  / volume(i,1));
        volumeDiffMin(i) = (min(volume(i,:))  / volume(i,1));
        volumeDiffMean(i)= (mean(volume(i,:))  / volume(i,1));
        
        % SLOPE (omitting first and last section
        slopeCells      = (rundata.store_etab(2:end-1,feed_starting_store:feed_ending_store-1)-circshift(rundata.store_etab(2:end-1,feed_starting_store:feed_ending_store-1),-1)) ./repmat(2*rundata.dx_array(2:end-1),1,feed_ending_store-feed_starting_store);
        slopeSpaceAvg        = sum(slopeCells(1:end-1,:))/(size(slopeCells,1)-1);
        slopeDiffEnd(i)	= (slopeSpaceAvg(end) / slopeSpaceAvg(1)) * 100;
        slopeDiffMax(i) = (max(slopeSpaceAvg)  / slopeSpaceAvg(1)) * 100;
        slopeDiffMin(i) = (min(slopeSpaceAvg)  / slopeSpaceAvg(1)) * 100;
        slopeDiffMean(i)= (mean(slopeSpaceAvg)  / slopeSpaceAvg(1)) * 100;

        % GSD init
        %keyboard
        GSD_Sizes(i,:)       = rundata.gsd_initDsi;
        GSD_initSurfFsi(i,:) = rundata.gsd_initSurfFsi(1,:);          % GSD init surf:
        GSD_feedFsi(i,:)     = rundata.gsd_feedFsi(1,:);              % GSD Feed:
        
        %GSD Time series:
        % single values, take spatial mean of focusArray:
        GSD_TS_Dsg(i,:)  = mean(rundata.store_Dsg(focusArray,feed_starting_store:feed_ending_store-1));     % GSD Dsg:
        GSD_TS_90(i,:)   = mean(rundata.store_Ds90(focusArray,feed_starting_store:feed_ending_store-1));    % GSD D90:
        % distributions:
        GSD_TS_SurfFsi(i,:,:)  = mean(rundata.store_SurfFsi(focusArray,:,feed_starting_store:feed_ending_store-1));    % GSD subsurface (?):
        GSD_TS_TransPbi(i,:,:) = mean(rundata.store_TranspPbi(focusArray,:,feed_starting_store:feed_ending_store-1));    % GSD subsurface (?):
                
        if goThroughSubsurfaceLayers
            GSD_TS_Dssg(i,:) = mean(rundata.store_Dssg(focusArray,feed_starting_store:feed_ending_store-1));    % GSD subsurface (?):
        else
            disp('FAKE DSSG with using feed')
            nodes_N = length(rundata.nodes_loc);
            [Dssg,~,Dss90,~] = GSDComputing_vectorized_spatial(repmat(rundata.gsd_feedFsi(1,:),nodes_N,1),rundata.gsd_feedDsi,nGSD-1,nodes_N);
            GSD_TS_Dssg(i,:) = mean(repmat(Dssg(focusArray),1,feed_ending_store-feed_starting_store));    % GSD subsurface (?):
        end
        
        % FEED Volume
        %keyboard
        %feedVol_m         = sum((rundata.feed_TArray_qsi_const(feed_starting_step:feed_ending_step) + rundata.feed_TArray_qsi_pulsed(feed_starting_step:feed_ending_step))) * rundata.param_dt;
        %feedVol_m3_const(i)     = sum(rundata.feed_TArray_qsi_const(feed_starting_step:feed_ending_step)) * rundata.dx_array(1) * rundata.param_reachbankfull(1) * rundata.param_dt;
        %feedVol_m3_pulse(i)     = sum(rundata.feed_TArray_qsi_pulsed(feed_starting_step:feed_ending_step)) * rundata.dx_array(rundata.feednodes) * rundata.param_reachbankfull(rundata.feednodes) * rundata.param_dt;
        %feedVol_kg_const(i)     = feedVol_m3_const(i) * rundata.const_rhoS;
        %feedVol_kg_pulse(i)     = feedVol_m3_pulse(i) * rundata.const_rhoS;
        
        % Time when the next pulse would hit if the simulation time was
        % longer:
        TPulseDistance_store(i)    =  TPulseDelta_hArr(i) * 3600 / TstepDt / ToutputFreq;
        TPulseNextWouldHit_store(i) = floor(feed_starting_store_real + floor(NPulse(i)) * TPulseDistance_store(i));
        if TPulseNextWouldHit_store(i) > feed_ending_store
            TPulseNextWouldHit_store(i) = TPulseNextWouldHit_store(i) - TPulseDistance_store(i);
        end
        if RunConfig.MariaClaudiaPulses
            % just last 280 hours
            sizeOfUsedArray(i) =feed_ending_store - feed_starting_store;
        else
            % pulsed run
            sizeOfUsedArray(i) = TPulseNextWouldHit_store(i) - feed_starting_store;
        end
        % Qb Sediment Transport at 1 nodes from flume exit
        qb_out(i,:)       = rundata.store_qbx(focusExportNode,feed_starting_store:feed_ending_store-1);
        %qb_out_allnodes(i,:,:)       = rundata.store_qbx(:,feed_starting_store:feed_ending_store);
        qb_mean(i)          = mean(qb_out(i,1:sizeOfUsedArray(i)));
        qb_divcum(i,:)      = cumsum(qb_out(i,:) - qb_mean(i));
        %qb_divdiff(i,:)     = diff(qb_div(i,:));
        qb_output_m(i,:)    = cumsum(qb_out(i,:)) * rundata.param_dt;
        %qb_output_m3(i,:)   = cumsum(qb_out(i,:) .* rundata.dx_array(end) .* rundata.param_reachbankfull(end) * rundata.param_dt);
        %qb_output_kg(i,:)   = qb_output_m3(i,:) * rundata.const_rhoS;
        
        % Sediment Feed:
        qb_in_const_step = rundata.feed_TArray_qsi_const(feed_starting_step:feed_ending_step);
        qb_in_pulse_step = rundata.feed_TArray_qsi_pulsed(feed_starting_step:feed_ending_step);
        
        stepsInStore = length(qb_in_const_step)/length(qb_in_const_store);
        for j = 1:arrSize(2)
            a = (j-1) * stepsInStore + 1;
            b = a + stepsInStore - 1;
            qb_in_const_store(i,j) = mean(qb_in_const_step(a:b));
            qb_in_pulse_store(i,j) = mean(qb_in_pulse_step(a:b));
        end
        clear qb_in_const_step qb_in_pulse_step rundata.feed_TArray_qsi_const rundata.feed_TArray_qsi_pulsed
        
        % Average output per pulse-period
        
        nPulseBegins_store(i,1) = 1;
        pulse_timing_sec(i,1) = 1;
        n = 1;
        analysis_offset = (StartAnalysisHour * 3600 / TstepDt / ToutputFreq);
        while true
            nPulseEnds_store(i,n)     = round(nPulseBegins_store(i,n) + TPulseDistance_store(i));
            if (n > NPulse(i)) || (floor(nPulseEnds_store(i,n) - analysis_offset) > sizeOfUsedArray(i))
                break
            end
            % if we are in focus time frame:
            if nPulseBegins_store(i,n) > (StartAnalysisHour * 3600 / TstepDt / ToutputFreq)
                qb_mean_pulse(i,n)      = mean(qb_out(i,floor(nPulseBegins_store(i,n)-analysis_offset):floor(nPulseEnds_store(i,n)-analysis_offset)));
            end
            pulse_timing_sec(i,n+1) = pulse_timing_sec(i,n) + TPulseDelta_hArr(i) * 3600;
            nPulseBegins_store(i,n+1) = nPulseEnds_store(i,n);
            n = n +1;
        end
        lastPulse = n;
        %percentual QbDiffCum
        qb_divcumPerc(i,:) = (cumsum(qb_out(i,:) - qb_mean(i))*100)./qb_mean(i);
        %qb_input(i,:)     = rundata.feed_TArray_qsi(feed_starting_step:feed_ending_step);
        %NPrint = rundata.NPrint
        
        % slope as mean of spatial nodes:
        slope(i,:)      = mean(rundata.store_slope(focusArray,feed_starting_store:feed_ending_store-1));
        % slope as mean of time (over whole run)
        
        PBegins_anastore = min(nPulseBegins_store(i,lastPulse)-analysis_offset,size(slope,2));
        PEnds_anastore   = min(nPulseEnds_store(i,lastPulse)-analysis_offset-1,size(slope,2));
        
        if size(slope,2) < nPulseEnds_store(i,lastPulse)-analysis_offset-1
            error('!!!! Pulse Period does not end on last index !!!!s')
        end
        
        slope_mean(i)   = mean(slope(i,PBegins_anastore:PEnds_anastore));

        % ustar as spatial mean:
        ustar(i,:)      = mean(rundata.store_ustar(focusArray,feed_starting_store:feed_ending_store-1));
        % ustar as temporal mean (over whole run)
        ustar_mean(i)   = mean(ustar(i,PBegins_anastore:PEnds_anastore));
        
        const_rhoS = rundata.const_rhoS;
        param_reachbankfull = rundata.param_reachbankfull;
        dx_array = rundata.dx_array;
        
        % Drop run data:
        clear rundata
    end
end

% Min/Max values for plotting:
plotmax_qb_out      = max(max(qb_out(:,:)));
plotmin_qb_out      = min(min(qb_out(:,:)));
plotmax_qb_divcum   = max(max(qb_divcum(:,:)));
plotmin_qb_divcum   = min(min(qb_divcum(:,:)));
plotmax_qb_divcumPerc = max(max(qb_divcumPerc(:,:)));
plotmin_qb_divcumPerc = min(min(qb_divcumPerc(:,:)));
plotmax_volume      = max(max(volume(:,:)));
plotmin_volume      = min(min(volume(:,:)));
plotmax_qb_divdiff  = max(max(qb_divdiff(:,:)));
plotmin_qb_divdiff  = min(min(qb_divdiff(:,:)));
plotmax_slope       = max(max(slope(:,:)));
plotmin_slope       = min(min(slope(:,:)));
plotmax_ustar       = max(max(ustar(:,:)));
plotmin_ustar       = min(min(ustar(:,:)));

% Save workspace
save(WorkspaceFile);

clear all
close all
end

