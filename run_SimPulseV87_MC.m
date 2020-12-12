% Sample file that runs different pulse orders. Based on Paper 
% Müller, T., & Hassan, M. A. (2018). Fluvial response to changes in the magnitude and frequency of sediment supply in a 1-D model. Earth Surface Dynamics, 6(4), 1041-1057.
% also see
% Elgueta-Astaburuaga, M. A., & Hassan, M. A. (2019). Sediment storage, partial transport, and the evolution of an experimental gravel bed under changing sediment supply regimes. Geomorphology, 330, 1-12.

close all
clear all

startdir = pwd;
resultsDir  = fullfile(pwd,'Results');
outputDir   = resultsDir;
scratchDir  = resultsDir;

NumberOfRuns = 1;

addpath(startdir)
addpath('batchfiles')
addpath('modelfiles')
addpath(genpath('analysis'))

%GS_SigmaArray = [0.05 0.1 0.2 0.4 0.8 1.0 1.2 1.4 1.6 2.0]; %1.6 is close to M&C MBWC
GS_SigmaArray = [1.6]; %1.6 is close to M&C MBWC

% config differences for sensitivity calcs (shown in appendix of Tobias' PhD thesis)
SimConfig = {%
    '',   10, 2.0, 0.75, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.75, 0.50, 1.5, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.75, 0.55, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.75, 0.55, 1.5, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.75, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.75, 0.50, 1.5, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.75, 0.55, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.75, 0.55, 1.5, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.80, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.80, 0.50, 1.2, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.80, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.80, 0.50, 1.2, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.85, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.0, 0.85, 0.50, 1.2, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.85, 0.50, 1.0, char(['0';'C';'1';'4';'2';'C';'0']);...
    '',   10, 2.2, 0.85, 0.50, 1.2, char(['0';'C';'1';'4';'2';'C';'0'])};

for s=1:size(SimConfig,1)
    SimConfig(s,1) = {['V87_',num2str(SimConfig{s,2}),'s_t',...
        num2str(floor(SimConfig{s,3})),'p',num2str(floor(mod(SimConfig{s,3},1)*10)),...
        '_aup',num2str(SimConfig{s,4}*100),...
        '_aFp',num2str(SimConfig{s,5}*100),...
        '_nLa',num2str(floor(SimConfig{s,6})),'p',num2str(floor(mod(SimConfig{s,6},1)*10)),...
        '_OF_MCWC']}; %V85_1s_t1p0_aup65_aFp40_nLa2p0_OF_MCWC
end

% SimConfig = {'V82_1s_t1p0_aup75_aFp25_nLa1p0_OF_MCWC',   1, 1.0, 0.75, 0.25, 1.0, char(['0';'C';'1';'4';'2';'C';'0'])};
    

%% FOR MANUAL EXECUTION: STOP ON ERROR (NOT GOOD FOR WHOLE BATCH!!!)
%dbstop if error

RunConfig.outputStyle.dispSubtime   = true;

overwriteResults = false;
pickupruns = true;

% parfor starts a parallel loop. Change to normal 'for' to run one sim at a
% time.
parfor batchConf_i = 1:size(SimConfig,1)
    cd(startdir)
    [SimName] = genconf_VX_MCPulse(SimConfig{batchConf_i,1},...
        SimConfig{batchConf_i,2},...
        SimConfig{batchConf_i,3},...
        SimConfig{batchConf_i,4},...
        SimConfig{batchConf_i,5},...
        SimConfig{batchConf_i,6},...
        SimConfig{batchConf_i,7},...
        pickupruns);
    
    for sigma=1:length(GS_SigmaArray)
        %for sigma=length(GS_SigmaArray):-1:1
        cd(startdir)
        runConfigName = SimName{sigma,1};
        runConfigFile = strcat('runconfig_',runConfigName,'.mat');
        %mkdir(runConfigName)
        
        %% Execute simulations
        for i=1:NumberOfRuns
            runDir = strcat(resultsDir,'/',runConfigName,'/',runConfigName,'_',num2str(i));
            mkdir(runDir)
            copyfile(strcat(startdir,'/',runConfigFile),runDir);
            cd(runDir)
            runName = strcat(runConfigName,'_',num2str(i));
            if ~exist(strcat(runName,'.mat'),'file') || overwriteResults
                run_BESMo_wrapper(runConfigFile, i, scratchDir,overwriteResults);
            else
                disp(strcat('skipping simulation, as ',runName,'.mat exists'))
            end
        end
        %run_flumemodel_wrapper_analysis(runConfigFile,resultsDir, outputDir)
    end
end
