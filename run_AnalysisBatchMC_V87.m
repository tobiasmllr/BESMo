close all;
clear all;

startdir = pwd;

addpath('batchfiles')
addpath('modelfiles')
addpath(genpath('analysis'))

resultsDir  = fullfile(pwd,'Results');
outputDir   = resultsDir;

%GS_SigmaArray = [0.05 0.1 0.2 0.4 0.7 1.0 1.3 1.6 2.0]; %1.6 is close to M&C
GS_SigmaArray = [1.6]; %1.6 is close to M&C

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


TableName = 'V87_BatchTable.mat';

batch_time      = 280;  % in hours
batch_mass      = 2100; % in kg
batch_massIncrease  = 0;
batch_massCounts    = 1; % if > 1, then the batches are iterated with a mass increase of batch_massIncrease


% outputDir   = '~/Results_local';
% resultsDir  = '~/Results_local';

focus.run1=6;
focus.run2=8;
focus.run3=12;
focus.run4=15;

StartAnalysisHour = 0; % Start analysis in this hour of feed time!

% for debugging
dbstop on error
forceRemake = true;
forceRemakeBatch = true;
goThroughSubsurfaceLayers = false;

tableBatch = table();

for bname = 1:size(SimConfig,1)
    batchName_composite = strcat(char(SimConfig{bname,1}),int2str(batch_mass),'kg',int2str(batch_time),'h','Gsigp');
    
    %% Do single run analysis
    for sigma=1:length(GS_SigmaArray)
        %run_flumemodel_wrapper_analysis('runconfig_M9WC30000kg4000hGsigp300.txt')
        cd(startdir)
        runConfigName = strcat(batchName_composite,num2str(GS_SigmaArray(sigma)*1000));
        runConfigFile = strcat('runconfig_',runConfigName,'.mat');
        
        [tableOut] = run_flumemodel_wrapper_analysis(runConfigFile,resultsDir,outputDir,forceRemake, forceRemakeBatch, goThroughSubsurfaceLayers,StartAnalysisHour);
        tableOut.NameConfig = {runConfigName};
        tableOut.NameShort = {['\alpha=',num2str(SimConfig{bname,5}),', n_a=',num2str(SimConfig{bname,6})]};
    end
    tableBatch = [tableBatch; tableOut];

    %keyboard
    
    %tableBatch
    %% Do batch analysis
    % just use any runConfigFile, the main parameters are the same for all
    % different sigma and frequencies...
    %cd(startdir)
%     run_batch_wrapper_analysis_forMarwan(batchName_prefix,batch_mass,batch_time, batch_massIncrease, batch_massCounts, outputDir, GS_SigmaArray,forceRemake,focus);
end
cd(startdir)
save(TableName,'tableBatch');

close all
