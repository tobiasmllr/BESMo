function [tableOut] = run_BESMo_wrapper_analysis(nameconfigfile_batch, resultsDir, outputDir, forceRemake,forceRemakeBatch , goThroughSubsurfaceLayers,StartAnalysisHour)
% This file is executed by each worker individually

% Get input parameters from the input file:
%[configparams_batch] = read_inputfile(nameconfigfile_batch);

% Get input parameters from the input file:
load(nameconfigfile_batch);


% check if this is one of the old configurations
if exist('RunConfig')
    Fnames = fieldnames(RunConfig);
    for field=1:size(Fnames,1)
        RunParam.(Fnames{field}) = RunConfig.(Fnames{field});
    end
end

%RunParam = fkt_generateBatchConstants();
[RunParam] = fkt_generateBatch(RunParam);


%% ========================================
%  CREATE FOLDER STRUCTURE AND INPUT FILES:
% set up a folder name for the batch:
batchDirName = RunParam.NameOfBatch;
batchDirPath = strcat(resultsDir,'/',batchDirName);

% % if input number was interpreted as string, convert:
% if ~isnumeric(i)
%    i = str2num(i);
% end

for i = 1:length(RunParam.P_Number)
    current_run{1,i} = strcat(RunParam.NameOfBatch,'_',int2str(i));
end
% runDir = strcat(batchDirPath,'/',current_run);

% LOGFILE: Where to save the fprintf out?
%fdev = 1; % 1 gives terminal output, (does not work with parfor!)
% else: file output 'current_runi_log.txt'
fdev = fopen(strcat(RunParam.NameOfBatch,'.log'),'w');

% set up a folder name for the analysis:
analysisDirName = strcat(RunParam.NameOfBatch,'_analysis');
analysisDirPath = strcat(outputDir,'/',analysisDirName);

%% ========================================
%  DO BATCH ANALYSIS:

%cd(batchDirPath)

if forceRemakeBatch % prevent making new plots
    fprintf(fdev, 'Starting Batch Analysis...\n');
    disp(strcat('Analysing run: ',RunParam.NameOfBatch));
    [tableOut] = fkt_analyzeRuns(batchDirPath, current_run, analysisDirPath, RunParam, RunParam,fdev,forceRemake, goThroughSubsurfaceLayers,StartAnalysisHour);
else
    disp(strcat('Skipping analysing runs, no forceRemake: ',RunParam.NameOfBatch));
end

%% Zip result folders:
%cd(baseDir)
%fprintf(fdev, 'Zipping resulting folders...\n');
% zip batch run:
%zip(strcat(batchDirName,'.zip'),batchDirName);
% zip analysis:
%zip(strcat(analysisDirName,'.zip'),analysisDirName);

% Clear variables
clear RunParam
fclose(fdev);
% RESULT: is a .mat file in the run folder. I load these in an other
% script to access the variables saved in there.
end
