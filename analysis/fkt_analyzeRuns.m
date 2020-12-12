function [ tableOut ] = fkt_analyzeRuns(batchDir, current_run, analysisDir, RunConfig, RunParam, fdev, forceRemake, goThroughSubsurfaceLayers,StartAnalysisHour)
%FKT_ANALYZERUNS Summary of this function goes here
%   Detailed explanation goes here

format long g
tableOut = table();

%%
runcount = length(RunConfig.P_Number);

plotSingleRuns = 'subplotonly';
fig_visible = 'on';

PlotNameStart = RunConfig.NameOfBatch;
WorkspaceFile = strcat(PlotNameStart,'.mat');

if ~exist(analysisDir,'dir')
    mkdir(analysisDir)
end
cd(analysisDir)

%
disp(WorkspaceFile)
disp(analysisDir)

if ~exist(WorkspaceFile,'file') || forceRemake
    fprintf('Workspace file does not exist... calculating variables\n');
    fkt_analyzeRuns_importVar( current_run, batchDir, RunConfig, RunParam, WorkspaceFile, goThroughSubsurfaceLayers,StartAnalysisHour)
end

% load workspace
fprintf('Loading Workspace file ...\n');
load(WorkspaceFile);

%% PLOTTING SINGLE-RUN
TimeArraySecTotal = feed_starting_sec:RunParam.NtoPrint * RunParam.dt:feed_ending_sec;
TimeArrayHourTotal = TimeArraySecTotal ./ 3600;
%TimeArraySecFeed = 0:RunParam.NtoPrint * RunParam.dt:(feed_ending_sec-feed_starting_sec);
TimeArrayHourFeed = TimeArrayHourTotal;

% VirtVeloQuot variables
VirtVeloQuot = NaN(arrSize(1),RunConfig.P_Number(1));
pulse_timing_hour = NaN(arrSize(1),RunConfig.P_Number(1));
pulse_distance_hour = NaN(arrSize(1),RunConfig.P_Number(1));
qb_mean_pulseKGpS = NaN(arrSize(1),RunConfig.P_Number(1));
pulseMag_kg_npulse = NaN(arrSize(1),RunConfig.P_Number(1));
pulseVirtVelocityKGpS = NaN(arrSize(1),1);

% feed variables:
feed_totalKg = NaN(arrSize(1),RunConfig.P_Number(1));
feed_totalM3 = NaN(arrSize(1),RunConfig.P_Number(1));

% preparing structures:
Output = struct();
plotQbDivAll(runcount) = struct();
plotQbOutAll(runcount) = struct();
plotSlopeAll(runcount) = struct();

%% PRECALCULATIONS:
for i=1:runcount
    % VIRTUAL VELOCITIES
    pulseVirtVelocityKGpS(i) = RunConfig.P_magnitudeArr(i)/(RunConfig.P_distanceArrayInMin(i) * 60);
end

for i=1:runcount
    for n = 1:RunConfig.P_Number(i)
        % Calculations
        pulseVirtVelocityKGpS_npulse(i,n) = pulseVirtVelocityKGpS(i);
        pulseMag_kg_npulse(i,n)           = RunConfig.P_magnitudeArr(i);
        % VirtFluvi_mean_pulse_KGpS
        qb_mean_pulseKGpS(i,n)  = qb_mean_pulse(i,n) * (1-RunParam.Lambda) * RunParam.RhoS * RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode) * RunParam.dx_array(RunParam.pulsedfeednode);
        pulse_timing_hour(i,n)  = pulse_timing_sec(i,n) / 3600;
        pulse_distance_hour(i,n)= RunConfig.P_distanceArrayInMin(i) / 60;
        VirtVeloQuot(i,n)       = qb_mean_pulseKGpS(i,n)/pulseVirtVelocityKGpS_npulse(i,n);
        %feed_totalKg(i,n)       = feedVol_kg_pulse(i)+feedVol_kg_const(i);
        %feed_totalM3(i,n)       = feedVol_m3_pulse(i)+feedVol_m3_const(i);
    end
    
end

plot_fontsize = 8;
plot_position =  [0, 0, 10, 4]; % in inches
plot_paperposition = [0, 0, 10, 4]; % in inches
plot_papersize = [10, 4]; % in inches

% batch run loop
runs = 1:runcount; %BatchConfig.run_starting:BatchConfig.run_overstep:BatchConfig.run_ending;
for j=1:runcount
    i = j;
    runDir = strcat(batchDir,'/',current_run{1,i});
    filename_mat = strcat(batchDir,'/',current_run{1,i},'.mat');
    if ~exist(filename_mat,'file')
        filename_mat = strcat(runDir,'/',current_run{1,i},'.mat');
        if ~exist(filename_mat,'file')
            msg = strcat('Error: File not found: ',filename_mat);
            warning(msg)
        end
    end
    fprintf(strcat('Plotting of analysis mat-file: ',filename_mat,'... \n'))
    
    
    if exist(filename_mat,'file')
        %fprintf(fdev,strcat('Plotting run ',num2str(i),'... \n'));
        
        % plotting set up
        if ~strcmp(plotSingleRuns,'none')
            %keyboard
            %fig_visible = 'on';
            %%
            format long g
            fig = figure('visible',fig_visible);
            set(fig,'DefaultAxesFontSize',6);
            set(fig,'DefaultTextFontSize',6);
            
            set(fig,'Units','Inches');
            set(fig, ...
                'Position', [0, 0, 11, 8.5],...
                'PaperPosition', [0, 0, 11, 8.5], ...
                'PaperUnits', 'Inches', ...
                'PaperSize', [11, 8.5]);
            
            col = 6; % subplot
            row = 2; % subplot
            
            %% ############### QB_Out ###############
            % convert qb from m/s to g/m/s
            qb_out_mps  = qb_out .* RunParam.NtoPrint .* RunParam.dt;
            qb_out_m3ps = (qb_out_mps .* RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode) .* RunParam.dx_array(RunParam.pulsedfeednode));
            qb_out_gms = qb_out_m3ps .* RunParam.RhoS; % should be kg/m/s, not grams. somehow it's factor 1000 too low...
            
            if RunConfig.MariaClaudiaPulses
                DataOut = struct();
                DataOut.time = [1:1:size(qb_out,2)] .* RunParam.NtoPrint .* RunParam.dt ./ 3600;
                DataOut.qb_gms = squeeze(qb_out_gms(1,:));
                DataOut.Dg = squeeze(GSD_TS_Dsg(1,:));
                DataOut.D90 = squeeze(GSD_TS_90(1,:));
                DataOut.slope = squeeze(slope(1,:));
            end
            % CHECK IF VOLUME MAKES SENSE
            totalout_kg = nansum(qb_out_gms(1,:));
            
%             totalfeed_m  = sum(FeedArrayConst) + sum(FeedArrayPulsed);
%             totalfeed_m3 = (totalfeed_m .* RunParam.ChBankfullWidth_array(RunParam.pulsedfeednode) .* RunParam.dx_array(RunParam.pulsedfeednode));
%             totalfeed_kg = totalfeed_m3 .* RunParam.RhoS;
            %%
            plotQbAll = subplot(col,row,1:2);
            
            %figure
            TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            if RunConfig.MariaClaudiaPulses
                % split cumulative into 40h pulse phases
                DataArrayPlot = qb_out_gms(i,1:sizeOfUsedArray(i));
                plotXmin = TimeArrayPlot(1);
                plotXmax = TimeArrayPlot(end);
                plotYmin = min(DataArrayPlot);
                plotYmax = max(DataArrayPlot);
                
            else
                % pulsed run
                DataArrayPlot = qb_out_gms(i,1:sizeOfUsedArray(i));
                plotXmin = TimeArrayPlot(1);
                plotXmax = TimeArrayPlot(end);
                plotYmin = plotmin_qb_divcum;
                plotYmax = plotmax_qb_divcum;
            end
            
            PlotTitle = 'qb divcum';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'qb_{out} [g/m/s]';
            
            semilogy(TimeArrayPlot(:),DataArrayPlot(:),'-');
            xlabel(PlotXLabel)
            ylabel(PlotYLabel)
            title(PlotTitle)
            
            if RunConfig.MariaClaudiaPulses
                hold on;
                y_limits = [plotYmin plotYmax];
                % rectangle for feed runs:
                r1 = rectangle('position',[0 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r2 = rectangle('position',[80 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r3 = rectangle('position',[160 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r4 = rectangle('position',[240 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                uistack(r1, 'bottom')
                uistack(r2, 'bottom')
                uistack(r3, 'bottom')
                uistack(r4, 'bottom')
                set(gca,'layer','top')
            else
                hold off;
            end
            
            
            if RunConfig.MariaClaudiaPulses
                plotHandle = figure('Visible',fig_visible);
                fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                    plotXmin, plotXmax, plotYmin, plotYmax, ...
                    PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
                axis tight
                hold off;
                plotName   = strcat(PlotNameStart,'_run',int2str(i),'_TR');
                fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
                fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
            end
            
            %% ############### QB_DIVCUM ###############
            plotQbDivAll = subplot(col,row,3:4);
            %figure
            TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            if RunConfig.MariaClaudiaPulses
                % split cumulative into 40h pulse phases
                phaseLength = 40 * 3600 / RunParam.dt / RunParam.NtoPrint;
                phaseN = 7;
                phaseStart = 1;
                for phase=1:phaseN
                    %figure;
                    if phase > 1
                        phaseStart = phaseStart + phaseLength;
                    end
                    phaseEnd   = phaseStart + phaseLength - 1;
                    qb_divcum(i,phaseStart:phaseEnd) = cumsum(qb_out_gms(i,phaseStart:phaseEnd) - mean(qb_out_gms(i,phaseStart:phaseEnd)));
                    %plot(qb_divcum(i,phaseStart:phaseEnd))
                end
                DataArrayPlot = qb_divcum(i,1:sizeOfUsedArray(i));
                plotXmin = TimeArrayPlot(1);
                plotXmax = TimeArrayPlot(end);
                plotYmin = min(DataArrayPlot);
                plotYmax = max(DataArrayPlot);
                
            else
                % pulsed run
                DataArrayPlot = qb_divcum(i,1:sizeOfUsedArray(i));
                plotXmin = TimeArrayPlot(1);
                plotXmax = TimeArrayPlot(end);
                plotYmin = plotmin_qb_divcum;
                plotYmax = plotmax_qb_divcum;
            end
            
            PlotTitle = 'qb divcum';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'qb diversion (qb_{out} - qb_{mean}) [g/m/s]';
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            
            if RunConfig.MariaClaudiaPulses
                hold on;
                y_limits = [plotYmin plotYmax];
                % rectangle for feed runs:
                try
                r1 = rectangle('position',[0 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r2 = rectangle('position',[80 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r3 = rectangle('position',[160 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r4 = rectangle('position',[240 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                uistack(r1, 'bottom')
                uistack(r2, 'bottom')
                uistack(r3, 'bottom')
                uistack(r4, 'bottom')
                set(gca,'layer','top')
                catch err
                    warning('rect didnt work...')
                end
            else
                hold off;
            end
            
            %             % ############### qb_output_rate ###############
            %             plotQbOutAll = subplot(col,row,3);
            %             DataArrayPlot = qb_out(i,1:sizeOfUsedArray(i));
            %             TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            %             PlotTitle = 'qb out';
            %             PlotXLabel= 'Feeded Experiment Time (hours)';
            %             PlotYLabel= 'qb output rate (m/s)';
            %             plotXmin = TimeArrayPlot(1);
            %             plotXmax = TimeArrayPlot(end);
            %             plotYmin = plotmin_qb_out;
            %             plotYmax = plotmax_qb_out;
            %             fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
            %                 plotXmin, plotXmax, plotYmin, plotYmax, ...
            %                 PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            %             % FIT QB
            %             y_mean = ones(length(DataArrayPlot),1) * mean(DataArrayPlot);
            %             plot(TimeArrayPlot, y_mean, 'r-','LineWidth',1);
            %             hold off;
            %
            %% ############### Slope ###############
            
            plotSlopeAll = subplot(col,row,5:6);
            DataArrayPlot = slope(i,1:sizeOfUsedArray(i));
            TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            PlotTitle = 'slope';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'slope (m/m)';
            plotXmin = TimeArrayPlot(1);
            plotXmax = TimeArrayPlot(end);
            plotYmin = plotmin_slope;
            plotYmax = plotmax_slope;
            %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            %hold on;
            % FIT Slope
            %y_mean = ones(length(DataArrayPlot),1) * mean(DataArrayPlot);
            %plot(TimeArrayPlot, y_mean, 'r-','LineWidth',1);
            if RunConfig.MariaClaudiaPulses
                hold on;
                y_limits = [plotYmin plotYmax];
                % rectangle for feed runs:
                r1 = rectangle('position',[0 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r2 = rectangle('position',[80 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r3 = rectangle('position',[160 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r4 = rectangle('position',[240 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                uistack(r1, 'bottom')
                uistack(r2, 'bottom')
                uistack(r3, 'bottom')
                uistack(r4, 'bottom')
                set(gca,'layer','top')
            else
                hold off;
            end
            
            %             % ############### Stored Volume to transport Rate 2 ###############
            %             subplot(col,row,7);
            %             plotStepLength = ceil(TPulseDistance_store(i)/8);
            %             fkt_plotValuePairDots_singleRunValue(i, ... %i
            %                 TimeArrayHourFeed(1:plotStepLength:sizeOfUsedArray(i)), ...   %TimeArray
            %                 volume(i,1:plotStepLength:sizeOfUsedArray(i)), ...     %DataArrayX
            %                 qb_out(i,1:plotStepLength:sizeOfUsedArray(i)), ...   %DataArrayY
            %                 'volume vs. qb out', ...                %PlotTitle
            %                 'qb output rate (m/s)', ...             %PlotYLabel
            %                 'stored volume (m3)');                  %PlotXLabel
            %             hold off;
            %
            
            % ############### GRAIN DSG surface ###############
            plotGS_SURFACE_All = subplot(col,row,7:8);
            
            DataArrayPlot = GSD_TS_Dsg(i,1:sizeOfUsedArray(i));
            TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            PlotTitle = 'Grain Size surface';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'Grain Size surface geometric mean(mm)';
            plotXmin = TimeArrayPlot(1);
            plotXmax = TimeArrayPlot(end);
            plotYmin = min(DataArrayPlot(:));
            plotYmax = max(DataArrayPlot(:));
            %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            if RunConfig.MariaClaudiaPulses
                hold on;
                y_limits = [plotYmin plotYmax];
                % rectangle for feed runs:
                r1 = rectangle('position',[0 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r2 = rectangle('position',[80 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r3 = rectangle('position',[160 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r4 = rectangle('position',[240 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                uistack(r1, 'bottom')
                uistack(r2, 'bottom')
                uistack(r3, 'bottom')
                uistack(r4, 'bottom')
                set(gca,'layer','top')
            else
                hold off;
            end
            
            
            % ############### GRAIN D90 surface ###############
            plotGS_D90SURFACE_All = subplot(col,row,9:10);
            
            DataArrayPlot = GSD_TS_90(i,1:sizeOfUsedArray(i));
            TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            PlotTitle = 'Grain Size D90 surface';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'Grain Size D90 surface(mm)';
            plotXmin = TimeArrayPlot(1);
            plotXmax = TimeArrayPlot(end);
            plotYmin = min(DataArrayPlot(:));
            plotYmax = max(DataArrayPlot(:));
            %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            if RunConfig.MariaClaudiaPulses
                hold on;
                y_limits = [plotYmin plotYmax];
                % rectangle for feed runs:
                r1 = rectangle('position',[0 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r2 = rectangle('position',[80 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r3 = rectangle('position',[160 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                r4 = rectangle('position',[240 y_limits(1) 40 diff(y_limits)],'FaceColor',[0.5 .9 .9],'EdgeColor','none');
                uistack(r1, 'bottom')
                uistack(r2, 'bottom')
                uistack(r3, 'bottom')
                uistack(r4, 'bottom')
                set(gca,'layer','top')
            else
                hold off;
            end
            
            if RunConfig.MariaClaudiaPulses
                plotHandle = figure('Visible',fig_visible);
                fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                    plotXmin, plotXmax, plotYmin, plotYmax, ...
                    PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
                axis tight
                hold off;
                plotName   = strcat(PlotNameStart,'_run',int2str(i),'_D90Surf');
                fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
                fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
            end
            
            
            %% Last hours focus
            skipsteps = 1;
            if RunConfig.MariaClaudiaPulses
                shownLastNHours = 280;
            else
                shownLastNHours = 400;
            end
            startstep       = find(TimeArrayHourTotal==TimeArrayHourTotal(end)-shownLastNHours);
            TimeArrayPlot = TimeArrayHourTotal(startstep:skipsteps:end-1);
            
            % ############### GRAIN Dg focus ###############
            plotGS_DgSURFACE_Last = subplot(col,row,11);
            
            DataArrayPlot = GSD_TS_Dsg(i,startstep:skipsteps:end);
            %TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            PlotTitle = 'Grain Size Dg surface';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'Grain Size Dg surface(mm)';
            plotXmin = TimeArrayPlot(1);
            plotXmax = TimeArrayPlot(end);
            plotYmin = min(DataArrayPlot(:));
            plotYmax = max(DataArrayPlot(:));
            %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            axis tight
            hold off;
            
            if RunConfig.MariaClaudiaPulses
                plotHandle = figure('Visible',fig_visible);
                fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                    plotXmin, plotXmax, plotYmin, plotYmax, ...
                    PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
                axis tight
                hold off;
                plotName   = strcat(PlotNameStart,'_run',int2str(i),'_DgSurf');
                pause(2);
                fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
                fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
            end
            
            %
            %             % ############### GRAIN D90 focus ###############
            %             plotGS_D90SURFACE_Last = subplot(col,row,9);
            %
            %             DataArrayPlot = GSD_TS_90(i,startstep:skipsteps:end);
            %             %TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            %             PlotTitle = 'Grain Size D90 surface';
            %             PlotXLabel= 'Feeded Experiment Time (hours)';
            %             PlotYLabel= 'Grain Size D90 surface(mm)';
            %             plotXmin = TimeArrayPlot(1);
            %             plotXmax = TimeArrayPlot(end);
            %             plotYmin = min(DataArrayPlot(:));
            %             plotYmax = max(DataArrayPlot(:));
            %             %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %             %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            %             fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
            %                 plotXmin, plotXmax, plotYmin, plotYmax, ...
            %                 PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            %             hold off;
            
            
            % ############### slope focus ###############
            plotGS_Slope_Last = subplot(col,row,12);
            
            DataArrayPlot = slope(i,startstep:skipsteps:end);
            %TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            PlotTitle = 'slope';
            PlotXLabel= 'Feeded Experiment Time (hours)';
            PlotYLabel= 'slope (m/m)';
            plotXmin = TimeArrayPlot(1);
            plotXmax = TimeArrayPlot(end);
            plotYmin = min(DataArrayPlot(:));
            plotYmax = max(DataArrayPlot(:));
            %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                plotXmin, plotXmax, plotYmin, plotYmax, ...
                PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            axis tight
            hold off;
            
            if RunConfig.MariaClaudiaPulses
                plotHandle = figure('Visible',fig_visible);
                fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
                    plotXmin, plotXmax, plotYmin, plotYmax, ...
                    PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
                axis tight
                hold off;
                pause(2);
                plotName   = strcat(PlotNameStart,'_run',int2str(i),'_Slope');
                fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
                fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                    plotHandle, plotName);
            end
            %             % ############### GRAIN DSG subsurface###############
            %             plotGS_SubSURFACE_All = subplot(col,row,6);
            %
            %             DataArrayPlot = GSD_TS_Dssg(i,1:sizeOfUsedArray(i));
            %             TimeArrayPlot = TimeArrayHourFeed(1:sizeOfUsedArray(i));
            %             PlotTitle = 'Grain Size sub surface';
            %             PlotXLabel= 'Feeded Experiment Time (hours)';
            %             PlotYLabel= 'Grain Size subsurface geometric mean(mm)';
            %             plotXmin = TimeArrayPlot(1);
            %             plotXmax = TimeArrayPlot(end);
            %             plotYmin = min(DataArrayPlot(:));
            %             plotYmax = max(DataArrayPlot(:));
            %             %backgroundcolorArr = ones(length(DataArray)) * plotYmax;
            %             %backgroundcolorArr(qb_div(i,:)<0) = NaN;
            %             fkt_plotTimeSeries_singleRunValue(i, TimeArrayPlot, DataArrayPlot(:), ...
            %                 plotXmin, plotXmax, plotYmin, plotYmax, ...
            %                 PlotTitle, PlotYLabel, PlotXLabel,sizeOfUsedArray(i));
            %             hold off;
            
            
            % Set plot title
            figtitle = strrep(current_run{1,i}, '_', ' ');
%             subtitle(strcat(figtitle, ' Pulsedistance = ',num2str(RunConfig.P_distanceArrayInMin(i) / 60), ' h'), 'fontsize',14);
            
            %%
            %keyboard
            %% save figure
            runplotname_png = strcat(PlotNameStart,'_run',int2str(i),'.png');
            runplotname_pdf = strcat(PlotNameStart,'_run',int2str(i),'.pdf');
            print(runplotname_png,'-dpng');
            print(fig,runplotname_pdf,'-dpdf','-painters');
            pdflist(j) = cellstr(runplotname_pdf);
        end
        
        if strcmp(plotSingleRuns,'all')
            %% Plot single run figures
            % QB DIV
            subplotHandle = plotQbDivAll;
            subplotName   = strcat(PlotNameStart,'_run',int2str(i),'_QbDivAll');
            fkt_printSubplotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            fkt_printSubplotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            
            % plotQbOutAll
            subplotHandle = plotQbOutAll;
            subplotName   = strcat(PlotNameStart,'_run',int2str(i),'_QbOutAll');
            fkt_printSubplotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            fkt_printSubplotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            
            % plotSlopeAll
            subplotHandle = plotSlopeAll;
            subplotName   = strcat(PlotNameStart,'_run',int2str(i),'_SlopeAll');
            fkt_printSubplotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            fkt_printSubplotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                subplotHandle, subplotName);
            
            
            % ############### Stored Volume to transport Rate 3 ###############
            plotVolToTRAll = figure('visible',fig_visible);
            plotStepLength = ceil(TPulseDistance_store(i)/8);
            fkt_plotValuePairDots_singleRunValue(i, ... %i
                TimeArrayHourFeed(1:plotStepLength:sizeOfUsedArray(i)), ...   %TimeArray
                volume(i,1:plotStepLength:sizeOfUsedArray(i)), ...     %DataArrayX
                qb_out(i,1:plotStepLength:sizeOfUsedArray(i)), ...   %DataArrayY
                'volume vs. qb out', ...                %PlotTitle
                'qb output rate (m/s)', ...             %PlotYLabel
                'stored volume (m3)');                  %PlotXLabel
            hold off;
            
            % plotVolToTRAll
            plotHandle = plotVolToTRAll;
            plotName   = strcat(PlotNameStart,'_run',int2str(i),'_VolToTRAll');
            fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                plotHandle, plotName);
            fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                plotHandle, plotName);
            
            % ############### VirtVeloQuot as Timeseries Plot ###############
            plotVirtVeloQuotTS_i= figure('visible',fig_visible);
            fkt_plotTimeSeries_singleRunValue(i, ...
                pulse_timing_hour(i,:), ...         %TimeArray
                VirtVeloQuot(i,:), ...              %DataArrayX
                min(pulse_timing_hour(i,:)), ...    %minx
                max(pulse_timing_hour(i,:)), ...    %maxx
                min(VirtVeloQuot(i,:)), ...         %miny
                max(VirtVeloQuot(i,:)), ...         %maxy
                'virtual velocity quotient', ...    %PlotTitle
                'U_{fluv}/U_{pulse}', ...           %PlotYLabel
                'Experiment time (hours)', ...      %PlotXLabel
                sizeOfUsedArray(i));
            % plotVolToTR
            plotHandle = plotVirtVeloQuotTS_i;
            plotName   = strcat(PlotNameStart,'_run',int2str(i),'_VirtVeloQuotTS');
            fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                plotHandle, plotName);
            fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
                plotHandle, plotName);
        end
        
        %% Store values to Structure --> Table --> CSV
        %Output(j).NameOfBatch = NameOfBatch;                          % ExperimentBatchName
        Output(j).i = i;                                               % RunNumber
        Output(j).TotalMagnitudeKg = (RunConfig.P_magnitudeArr(i) * RunConfig.P_Number(i));     % TotalMagnitude
        Output(j).NumberOfPulses = RunConfig.P_Number(i);                          % NumberOfPulses
        Output(j).DistanceOfPulsesMin = RunConfig.P_distanceArrayInMin(i);         % DistanceOfPulses
        Output(j).MagnitudeEachKg = RunConfig.P_magnitudeArr(i);                    % MagnitudeEach
        
        % Transport virtual / modelled
        Output(j).PulseVirtualVelocityKgpS = pulseVirtVelocityKGpS(i); % PulseVirtualVelocity
        MeanQbAll       = mean(qb_out(i,1:sizeOfUsedArray(i)));
        
        %Output(j).VelocityQuotient = Output(j).QbFitKgpS / Output(j).PulseVirtualVelocityKgpS;              %Relation between measured Qb and pulse virtual velocity
        
        % Mean slope
        %Output(j).MeanSlopeAll = MeanSlopeAll;
        MeanSlopeAll    = mean(slope(i,1:sizeOfUsedArray(i)));
        
        % Mean qb
        %Output(j).MeanQbAllMpS = MeanQbAll;
        % Slope fit
        %Output(j).SlopefitAllA = PolyFitSlopeAll(1);
        %Output(j).SlopefitAllB = PolyFitSlopeAll(2);
        %Output(j).SlopefitAllC = PolyFitSlopeAll(3);
        
        %Qb Fit
        %Output(j).QbfitAllA = PolyFitQbAll(1);
        %Output(j).QbfitAllB = PolyFitQbAll(2);
        %Output(j).QbfitAllC = PolyFitQbAll(3);
    end
end

try
    %% SAVE CSV TABLE
    filename = strcat(RunConfig.NameOfBatch,'_all.csv');
    OutputTable = struct2table(Output);
    writetable(OutputTable,filename,'WriteRowNames',true,'Delimiter',';');
    
    tableOut = table;
    tableOut.D90 = GSD_TS_90(1,:);
    tableOut.Dg  = GSD_TS_Dsg(1,:);
    tableOut.Slope  = slope(1,:);
    tableOut.TR_gms  = qb_out_gms(1,:);
    tableOut.time  = TimeArrayPlot(1,:);
catch err
    warning('Some value not found!!')
end
%% PLOT MULTI RUN

%keyboard
%save example workspace:
%save(strcat(PlotNameStart,'.mat'));
%
% %% ############### Virtual Velocity Quotient vs Mass ###############
% %fig_visible = 'off';
% figVirtVeloQuotvsMass = figure('visible',fig_visible);
%
% % Plotting
% fkt_plotValuePairDots_multiRunValue(...
%     pulse_timing_hour, ...                      %TimeArray
%     pulseMag_kg_npulse, ...                     %DataArrayX
%     VirtVeloQuot, ...                           %DataArrayY
%     'virtual velocity quotient', ...            %PlotTitle
%     'U_{fluv}/U_{pulse}', ...                   %PlotYLabel
%     'Pulse magnitude M_{pulse} [kg/pulse]',...  %PlotXLabel
%     'Experiment time (hours)',...               %PlotZLabel
%     runs,...                % run config
%     RunConfig.P_Number,...              % Number of Pulses
%     'true');                % Plot 1-line, 0.9-line and 0.8-line?
%
% % export to EPS and TEX
% plotHandle = figVirtVeloQuotvsMass;
% plotName   = strcat(PlotNameStart,'_all_VirtVeloQuotvsMass');
% fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
% fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
%
% %% ############### Virtual Velocity Quotient vs Pulse Distance ###############
% %fig_visible = 'on';
% figVirtVeloQuotvsPDist = figure('visible',fig_visible);
%
% % Plotting
% fkt_plotValuePairDots_multiRunValue(...
%     pulse_timing_hour, ...                      %TimeArray
%     pulse_distance_hour, ...                    %DataArrayX
%     VirtVeloQuot, ...                           %DataArrayY
%     'virtual velocity quotient', ...            %PlotTitle
%     'U_{fluv}/U_{pulse}', ...                   %PlotYLabel
%     'Pulse distance T_{pulse} [hours]',...      %PlotXLabel
%     'Experiment time (hours)',...               %PlotZLabel
%     runs,...                % run config
%     RunConfig.P_Number,...              % Number of Pulses
%     'true');                % Plot 1-line, 0.9-line and 0.8-line?
%
% % export to EPS and TEX
% plotHandle = figVirtVeloQuotvsPDist;
% plotName   = strcat(PlotNameStart,'_all_VirtVeloQuotvsPDist');
% fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
% fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
%
% %% ############### Virtual Velocity Quotient as Timeseries ###############
% %fig_visible = 'on';
% plotVirtVeloQuotTS= figure('visible',fig_visible);
%
% % Plotting
% fkt_plotTimeSeries_multiRunValue( ...
%     pulse_timing_hour, ...              %TimeArray
%     VirtVeloQuot, ...                   %DataArrayX
%     pulseMag_kg_npulse(:,1), ...        %LegendItems_iArray1
%     'kg', ...                           %LegendItemText1
%     pulse_distance_hour(:,1), ...       %LegendItems_iArray2 (will be used if all items of Legend 1 are the same)
%     'h', ...                            %LegendItemText2
%     'virtual velocity quotient', ...    %PlotTitle
%     'U_{fluv}/U_{pulse}', ...           %PlotYLabel
%     'Experiment time (hours)',...       %PlotXLabel
%     'true', ...                         %Plot 1-line, 0.9-line and 0.8-line?
%     20);                                %Maximum number of legend entries if more lines, they will be spread
%
% % export to EPS and TEX
% plotHandle = plotVirtVeloQuotTS;
% plotName   = strcat(PlotNameStart,'_all_VirtVeloQuotTS');
% fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
% fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
%
% %% Initial GSD
% %fig_visible = 'off';
%
% plotInitGSD= figure('visible',fig_visible);
% plot(GSD_initSurfFsi(1,:))
%
% matlabver = version('-release');
% if ~strcmp(matlabver,'2013b')
%     set(gca,'XTickLabelRotation',45);
% end
%
% xlabel('Grain size (mm)') % x-axis label
% ylabel('Frequency') % y-axis label
% set(gca,'XTick',[1:17],...
%     'XTickLabel',GSD_Sizes(1,:));
%
% hold off;
%
% % export to EPS and TEX
% plotHandle = plotInitGSD;
% plotName   = strcat(PlotNameStart,'_all_initialGSD');
% fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
% fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)

% %% ############### GSD as Timeseries ###############
% fig_visible = 'on';
% plotGSDTS= figure('visible',fig_visible);
%
% % Plotting
% fkt_plotTimeSeries_multiRunValue( ...
%     pulse_timing_hour, ...              %TimeArray
%     GSD_TS_Dsg, ...                   %DataArrayX
%     pulseMag_kg_npulse(:,1), ...        %LegendItems_iArray1
%     'kg', ...                           %LegendItemText1
%     pulse_distance_hour(:,1), ...       %LegendItems_iArray2 (will be used if all items of Legend 1 are the same)
%     'h', ...                            %LegendItemText2
%     'virtual velocity quotient', ...    %PlotTitle
%     'U_{fluv}/U_{pulse}', ...           %PlotYLabel
%     'Experiment time (hours)',...       %PlotXLabel
%     'true', ...                         %Plot 1-line, 0.9-line and 0.8-line?
%     20);                                %Maximum number of legend entries if more lines, they will be spread
%
% %% export to EPS and TEX
% plotHandle = plotGSDTS;
% plotName   = strcat(PlotNameStart,'_all_GSDTS');
% fkt_printPlotToEPS( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)
% fkt_printPlotToTEX( fig_visible, plot_fontsize, plot_position, plot_paperposition, plot_papersize, ...
%     plotHandle, plotName)

%% SAVE ALL RUNS PDF (error?):
if ~strcmp(plotSingleRuns,'none')
    batchplotname = strcat(RunConfig.NameOfBatch,'_all.pdf');
    if exist(batchplotname,'file')
        delete(batchplotname);
    end
    if exist('pdflist')
        try
            append_pdfs(batchplotname, pdflist{:});
        catch
            warning('might need https://www.mathworks.com/matlabcentral/fileexchange/31215-append_pdfs')
        end
    end
end

% % QB_DIFFCUMPERC
% DataArray = qb_diffcumPerc(:,:);
% PlotTitle = 'qb_diffcumPerc';
% PlotXLabel= 'Feeded Experiment Time (hours)';
% PlotYLabel= 'qb difference (qb_out - qb_mean)*100/qb_mean';
% plotYmin   = plotmin_qb_diffcumPerc;
% plotYmax   = plotmax_qb_diffcumPerc;
% %fig_visible = 'off';
% for i=run_starting:run_overstep:run_ending
%     current_run = strcat(NameOfBatch,num2str(i));
%     filename = strcat(offlinedata_dir,current_run,'.mat');
%     if exist(filename)
%         fkt_plotTimeSeries_singleRunValue(i, plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%             TimeArrayHour, DataArray(i,:), plotXmin, plotXmax, plotYmin, plotYmax, PlotTitle, PlotYLabel, PlotXLabel,PlotNameStart,sizeOfUsedArray(i))
%     end
% end

% QB_DIFFCUMPERC_ALLInOne
%DataArray = qb_divcumPerc(:,:);
%PlotTitle = 'qb_diffcumPerc';
%PlotXLabel= 'Feeded Experiment Time (hours)';
%PlotYLabel= 'qb difference (qb_out - qb_mean)*100/qb_mean';
%fig_visible = 'off';
%fkt_plotTimeSeries_singleRunValueOnePlot(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%    TimeArrayHour, DataArray, PlotTitle, PlotYLabel, PlotXLabel,PlotNameStart,NameOfBatch,offlinedata_dir,sizeOfUsedArray(i))

% % qb_output_rate ALL NODES
% DataArray = qb_out_allnodes(:,:,:);
% PlotTitle = 'qb_out_allnodes';
% PlotXLabel= 'Feeded Experiment Time (hours)';
% PlotYLabel= 'qb output rate all nodes';
% for i=1:size(DataArray,1)
%     current_run = strcat(NameOfBatch,num2str(i));
%     filename = strcat(offlinedata_dir,current_run,'.mat');
%     if exist(filename)
%         fkt_plotTimeSeries_nodesRunValue(i, plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%             TimeArrayHour, TstepDt, ToutputFreq, DataArray, PlotTitle, PlotYLabel, PlotXLabel,PlotNameStart)
%     end
% end


% % qb_output_kg
% TimeArraySec = 1:TstepDt*ToutputFreq:(size(qb_out,2))*TstepDt*ToutputFreq;
% TimeArrayHour = TimeArraySec ./ 3600;
% DataArray = qb_output_kg;
% PlotTitle = 'qb_output_kg';
% PlotXLabel= 'Feeded Experiment Time (hours)';
% PlotYLabel= 'qb output (kg)';
% for i=1:size(DataArray,1)
%     current_run = strcat(NameOfBatch,num2str(i));
%     filename = strcat(offlinedata_dir,current_run,'.mat');
%     if exist(filename)
%         fkt_plotTimeSeries_singleRunValue(i, plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%             TimeArrayHour, DataArray(i,:), PlotTitle, PlotYLabel, PlotXLabel,PlotNameStart)
%     end
% end

%% PLOTTING MULTI-RUN
% TimeArrayHour = TPulseDelta_hArr;
%
% %SLOPE
% DataArrayEnd    = slopeDiffEnd;
% DataArrayMin    = slopeDiffMin;
% DataArrayMean   = slopeDiffMean;
% DataArrayMax    = slopeDiffMax;
% PlotTitle       = 'Slope';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'Slope difference (% of start)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
% %D90 in Relation to Bed initial
% DataArrayEnd    = GSD90End;
% DataArrayMin    = GSD90Min;
% DataArrayMean   = GSD90Mean;
% DataArrayMax    = GSD90Max;
% PlotTitle       = 'D90';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'D90 (% of start)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
%
% %D90 in Relation to Feed
% DataArrayEnd    = GSD90EndToFeed;
% DataArrayMin    = GSD90MinToFeed;
% DataArrayMean   = GSD90MeanToFeed;
% DataArrayMax    = GSD90MaxToFeed;
% PlotTitle       = 'D90tofeed';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'D90 (% of feed)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
%
% %DSG in relation to initial bed
% DataArrayEnd    = GSDSGEnd;
% DataArrayMin    = GSDSGMin;
% DataArrayMean   = GSDSGMean;
% DataArrayMax    = GSDSGMax;
% PlotTitle       = 'Dsg';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'Dsg (% of start)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
% %DSG in relation to feed
% DataArrayEnd    = GSDSGEndToFeed;
% DataArrayMin    = GSDSGMinToFeed;
% DataArrayMean   = GSDSGMeanToFeed;
% DataArrayMax    = GSDSGMaxToFeed;
% PlotTitle       = 'Dsgtofeed';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'Dsg (% of feed)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
% %VOLUME Difference
% DataArrayEnd    = volumeDiffEnd;
% DataArrayMin    = volumeDiffMin;
% DataArrayMean   = volumeDiffMean;
% DataArrayMax    = volumeDiffMax;
% PlotTitle       = 'Volume';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'Volume (% of start)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
%
% % Feed Mass Total
% DataArrayEnd    = feedVol_kg;
% DataArrayMin    = NaN;
% DataArrayMean   = NaN;
% DataArrayMax    = NaN;
% PlotTitle       = 'FeedMassTotal';
% PlotXLabel      = 'Pulse distance (hours)';
% PlotYLabel      = 'Mass (kg)';
% fkt_plotTimeSeries_EndMinMeanMax(plot_x, plot_y, plot_width, plot_height, run_index_tru, fig_visible, ...
%     TimeArrayHour, DataArrayEnd, DataArrayMin, DataArrayMean, DataArrayMax, PlotTitle, PlotYLabel, PlotXLabel, PlotNameStart)
end

