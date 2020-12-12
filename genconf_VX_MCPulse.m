function [NameOfBatch] = genconf_VX_MCPulse(batchName_prefix,timestep,tau_crit_factor,param_au,alfaF,nLa_activeLayer,pulseOrder,pickupruns)
% GENERATE MULTI BATCH RUN PARAMETERS:
% 1: Type
%   A: Total Magnitude devided over N Pulses
%   B: Pulses of constant magnitude over N Pulses
%   C: Pulses of constant magnitude with different temporal pulse distance
%
% 2: Magnitude(kg)
%   if Type A: Total Magnitude
%   if Type B: Magnitude per Pulse
%   if Type C:
%
% 3: Time without feed before experiment starts (in hours):
% 4: Feed duration (in hours)
% 5: Time without feed after experiment ends (in hours):
%   3 + 4 + 5 = total simulation time
%
% 6: Pulse distance (start) (in minutes)
% 7: Pulse distance (increase) (in minutes)
% 8: Pulse duration (time it takes to feed the pulse mag in each pulse) (in minutes)
% 9: Pulse N minimum (N over feed time)
%
%10: GSD size classes index :
%   if 1: [0.25 1 4]        2 Grain size classes (in mm)
%   if 2: [0.25 1 2 4 8]    4 Grain size classes (in mm)
%   if 3: [0.177 0.25 0.354 0.5 0.71 1 1.41 2 2.83 4 5.66 8 11.2 16 22.6 32 44] (in mm)
%   but see fkt_generateBatch to see exact indices!
%11: GSD percent finer (initial bed):
%   if 1: [0 0 100]
%   if 2: [0 20 40 80 100]
%   if 3: [0 1.971 3.942 9.363 15.276 22.681 30.566 41.407 51.262 61.117 71.410 81.702 89.136 93.138 95.426 98.856 100.000]
%12: GSD percent finer (sediment feed):
%   see option 9...
%13: Water discharge (constant over run, in m3/s)
%14: Initial bed slope (in m/m)
%15: BatchName
%% Set variables
format long e
addpath('batchfiles')
addpath('modelfiles')

RunParam.slopeOrder = 1;
RunConfig.continuePreviousConfig = false;

%% SIMULATION TIME AND FEED PARAMETERS:
RunConfig.T_before = 0; 
RunConfig.T_sim = 280; %18000 ~2 years
RunConfig.T_after = 0;
RunConfig.T_total = RunConfig.T_before + RunConfig.T_sim + RunConfig.T_after; % in hours

RunConfig.P_distanceStart = 10; 
RunConfig.P_distanceIncrease = 20; % for A and B: increase of pulse distance between pulses in min. For Log: Number of Pulses simulated (log10 distributed in distances)
RunConfig.P_duration = 10; %pulse duration in minutes
PNMaxHourPerPulse = 150; %Max Pulse distance in hours!
RunConfig.P_NMin = floor(RunConfig.T_sim/PNMaxHourPerPulse); % Minimum Number of Pulses per Run

RunConfig.RunType = 'An'; 
% batchType A: Split simulation time into pulses, change magnitude (same mass fed/time)
% batchType An: Same as A, but only keep PDIncr number of pulses that are more evenly distributed between PDInit und PNMaxHourPerPulse;
% batchType B: Just change pulse distance within batch, constant magnitude (different mass fed/time)

batchCount = 1;

FeedMass = 0.5 * RunConfig.T_sim * (300/40);
FeedMassIncreasePerBatch = 0.5 * RunConfig.T_sim * (300/40);

% Recreate Maria and Claudia's experiment pulse distribution, ignoring some
% of the parameters here and overwriting them in the pulse generation
% script later on...
RunConfig.MariaClaudiaPulses = true;
%RunConfig.MariaClaudiaSequence = char(['0';'C';'1';'4';'2';'C';'0']);
%RunConfig.MariaClaudiaSequence = char(['0';'C';'4';'2';'1';'C';'0']);
%RunConfig.MariaClaudiaSequence = char(['0';'C';'1';'2';'4';'C';'0']);
%RunConfig.MariaClaudiaSequence = char(['C';'4';'C';'2';'C';'1';'C']);
%RunConfig.MariaClaudiaSequence = char(['C';'1';'C';'2';'C';'4';'C']);
RunConfig.MariaClaudiaSequence = pulseOrder;

RunParam.constantFeedFract = 0; % fraction of total feed that is fed constantly over sim time
RunParam.tau_crit_factor = tau_crit_factor;

filename_geometry = strcat('geometry_flume1_50cm.csv');
%% SEDIMENT TRANSPORT EQUATION
RunParam.eqNumber = 2;
%          1: Parker (1990)
%          2: Wilcock-Crowe (2003)
%          3: Ashida-Michiue (1972)
%          4: Meyer-Peter Müller

%% GRAIN SIZE PARAMETERS
% GS_Classes =  1; %: D50
% GS_Classes =  3; %: D16 D50 D84
% GS_Classes =  5; %: D10 D16 D50 D84 D90
% GS_Classes =  7; %: D10 D16 D34 D50 D66 D84 D90
% GS_Classes =  9; %: D5 D10 D16 D34 D50 D66 D84 D90 D95
RunConfig.GS_Classes = 11; %: D5 D10 D16 D25 D34 GS_D50 D66 D75 D84 D90 D95

% Calculate grain size distribution
RunConfig.GS_D50 = 5.640394298438923e+00;

%GS_SigmaArray = [0.05 0.1 0.2 0.4 0.7 1.0 1.3 1.6 2.0]; %1.6 is close to M&C
GS_SigmaArray = [1.6]; %1.6 is close to M&C

RunConfig.plotStyle = {'b-','b--','m-','m--','r-','r--','k-','k--','c-','c--','g-','g--','y-','y--'};

RunConfig.GSDistributionType = 'Normal'; % alternative: 'LogNormal'

%% DISCHARGE PARAMETERS
RunConfig.Discharge      = 0.065; %m2/s
RunParam.freqQ          = 1;    % freqQ: Hydrograph frequency

%% TEMPORAL PARAMETERS
RunParam.dt             = timestep;    % dt: time step, in s
RunParam.NtoPrint       = 10 * 60 / RunParam.dt;   % NtoPrint: Number of steps until a printout is made

RunParam.NPrint = RunConfig.T_total * 3600 / RunParam.NtoPrint / RunParam.dt;  

%% PHYSICAL CONSTANTS
%PhysicalConstants.txt
RunParam.RhoS           = 2650;     % rhos: Sediment density, in kg/m3
RunParam.Lambda         = 0.4;      % lambda: bed porosity
RunParam.Gravity        = 9.81;     % g: Acceleration of gravity, in m/s2
RunParam.RhoW           = 1000;     % rho: Water density, in kg/m3

%% AUXILIARY PARAMETERS
%AuxiliaryParameters.txt
RunParam.nk             = 2;        % nk: Dimensionless Nikuradse number
RunParam.alfar          = 8.1;      % alfar: Dimensionless Manning - Strickler coefficient evaluated by means of Parker (1991)
RunParam.alfaF          = alfaF;      % alfaF: Dimensionless coefficient for exchange fractions between the Active Layer and the Substrate
RunParam.nLa            = nLa_activeLayer;        % nLa: Dimensionless coefficient for the active layer thickness

RunParam.useUpwinding   = true;     % CAUSES PROBLEMS WITH FEEDING IN THE MIDDLE OF THE SIM SPACE
RunParam.au             = param_au;     % au: Upwinding coefficient


%StratigraphyData.txt
RunParam.deta_s         = 0.1;   % deta_s = Thickness of the layer to store the stratigraphy, in m 
% alfa_LZu changed from 10 to 25 to accomodate large pulses: /tobias
RunParam.alfa_LZu       = 30;     % alfa_LZu =  Dimensionless parameter that define top bed deposit elevation to consider in the stratigraphy
RunParam.alfa_LZd       = 10.0;    % alfa_LZd =  Dimensionless parameter that define bottom bed deposit elevation to consider in the stratigraphy

% parameters for backwater function
RunParam.backwaterparam.iter_max       = 1000;      % maximum number of iterations to reach iter_diffgood:
RunParam.backwaterparam.iter_diffgood  = 0.0001;    % stop iterating waterdepth approximation if difference between calculations is smaller than this value
RunParam.backwaterparam.waterdepth_min = 0.0001;    % calc finer spatial grid if waterdepth is calculated to be smaller than value. maybe calculate this threshold realistically?
RunParam.backwaterparam.froudefact_min   = 0.5;     % calc finer spatial grid if upstream decrease of froudenumber in subcritical condition is calculated to be smaller than this value multiplied by downstream Froude number.
RunParam.backwaterparam.subGridScalemax = 65536;    % maximum interpolation factor for calculating subgrids. Error if tripped
RunParam.backwaterparam.Fr_threshold     = 0.9;     % Froude Number Threshold for using supercritical:


%% SET NETCDF PARAMETERS
RunConfig.nc_params = struct();
% Using netcdf:
RunConfig.nc_params.saveNC        = true;
% saving subsurface data greatly increases processing time! I use it only
% for video output... Should be used for being able to pick up runs,
% though! Consider increasing NtoPrint to store less frequently
RunConfig.nc_params.saveNC_subs   = true;

% parameters
RunConfig.nc_params.deflate       = true;
RunConfig.nc_params.shufflefilter = false;
RunConfig.nc_params.setFill       = 'FILL';
RunConfig.nc_params.syncInterval  = 50;     % sync data to disk every index * store_netcdf_interval.

% Chunk size parameters:
%   I found these parameters through testing... They greatly affect storage
%   performance
%   I got crashes with premp = 1;
RunConfig.nc_params.deflate_level = 1;
RunConfig.nc_params.csize     = 45000000; % in bytes
RunConfig.nc_params.nelems    = 1500; % number of elements
RunConfig.nc_params.premp     = 0.1; % fraction preempt (between 0 and 1)

%% SET OUTPUT PARAMETERS-
RunConfig.outputStyle.dispSubtime   = true;
RunConfig.outputStyle.showStoreSize = true;

%% SET RUN MODE // DEBUG OR 
RunConfig.run_mode.debug              = false;  % if true         , run exits after index==debugSteps
RunConfig.run_mode.debugSteps         = 500;    % if debug is true, run exits after index==debugSteps
RunConfig.run_mode.debugGenFilename   = false;  % if true, generate filename for nc from individual time stamp, instead of current_run

% try to prevent negative slopes with introducing subtimesteps
RunConfig.run_mode.doSubtimesteps     = false;

% try to pick up old run from nc? 
%NOT IMPLEMENTED RIGHT NOW
RunConfig.run_mode.pickup             = pickupruns;

%% GEOMETRY PARAMETERS (NEW STUFF FOR P1WC and on, input geometry from csv)
% Use input geometry from csv OR if false:
% calculate node positions and slope/elevation assuming straight starting 
% geometry (might have problem with backwater flow causing negative slope...
RunParam.useGeometryCSV =false;

if ~RunParam.useGeometryCSV
    % possible problem with negative slope due to backwater flow
    dx              = 1;  % dx: cell size, in m
    Slope0          = 0.022;
    chwidth           = 1;
    
    RunParam.ReachLength     = 12;        % L: Length of the reach, in m
    RunParam.nodes_loc       = 0 : dx : RunParam.ReachLength;
    RunParam.ChBankfullWidth_array = ones(1,RunParam.ReachLength+1) .* chwidth;        % Channel Bankfull width, in m
    RunParam.slope_array     = ones(length(RunParam.nodes_loc),1) .* Slope0;
    RunParam.dx_array        = ones(length(RunParam.nodes_loc),1) .* dx;
    RunParam.elev_array      = RunParam.slope_array .* transpose(RunParam.ReachLength - RunParam.nodes_loc);
    
    RunParam.constfeednode = 1;
    RunParam.pulsedfeednode = 1;
else
    % read geometry input file:
    inputtable = readtable(filename_geometry);

    % get geometry
    RunParam.nodes_loc              = table2array(inputtable(:,2));
    RunParam.dx_array               = NaN(length(RunParam.nodes_loc),1);
    RunParam.dx_array(1:end-1)      = RunParam.nodes_loc(2:end) - RunParam.nodes_loc(1:end-1);
    RunParam.dx_array(end)          = RunParam.dx_array(end -1);
    RunParam.ReachLength            = RunParam.nodes_loc(end)-RunParam.nodes_loc(1);
    RunParam.elev_array             = table2array(inputtable(:,3));
    RunParam.slope_array            = RunParam.elev_array./(RunParam.ReachLength-RunParam.nodes_loc);
    RunParam.slope_array(end)       = RunParam.slope_array(end-1); % get slope for last node
    RunParam.ChBankfullWidth_array  = table2array(inputtable(:,4));        % Channel Bankfull width, in m
    RunParam.erodability_array      = table2array(inputtable(:,5));
    
    % get feed positions
    constfeednode_array     = table2array(inputtable(:,6));
    pulsedfeednode_array     = table2array(inputtable(:,7));
    RunParam.constfeednode = find(constfeednode_array,1,'first');
    RunParam.pulsedfeednode = find(pulsedfeednode_array,1,'first');
    
    % check input parameters:
    assert(sum(constfeednode_array)==1,             'Number of defined constant feed nodes is not 1.');
    assert(sum(pulsedfeednode_array)==1,            'Number of defined pulsed feed nodes is not 1.');
    assert(~any(RunParam.ChBankfullWidth_array<=0), 'There are nodes with width <= 0');
    assert(~any(RunParam.elev_array<0),             'There are nodes with elevation < 0');
    assert(~any(RunParam.dx_array<=0),              'There are nodes with spacing <= 0');
end

%% Generate run batch name:
for i=1:length(GS_SigmaArray)
    batchName_suffix = strcat('kg',num2str(RunConfig.T_sim),'h','Gsigp',num2str(GS_SigmaArray(i)*1000));

    for j=1:batchCount
        mass = FeedMass + (j * FeedMassIncreasePerBatch);
        if mass < 10
            NameOfBatch{i,j} = strcat(batchName_prefix,'000', num2str(mass), batchName_suffix);
        elseif mass < 100
            NameOfBatch{i,j} = strcat(batchName_prefix,'00', num2str(mass), batchName_suffix);
        elseif mass < 1000
            NameOfBatch{i,j} = strcat(batchName_prefix,'0', num2str(mass), batchName_suffix);
        else
            NameOfBatch{i,j} = strcat(batchName_prefix,num2str(mass), batchName_suffix);
        end
        
        % If type A: total mass fed, If type B: mass per pulse
        %batchType A: total mass fed in grams [300kg/40h ^= 1500kg/200h ^= 3000kg/400h ^= 15000kg/2000h ^= 30000]
        %batchType B: mass per pulse in grams
        RunConfig.P_magnitude       = mass;
        
        RunConfig.GS_Sigma          = GS_SigmaArray(i);
        
        RunConfig.NameOfBatch = NameOfBatch{i,j};
        %% CREATE RUNCONFIG mat
        filename = strcat('runconfig_',NameOfBatch{i,j},'.mat');
        save(filename,'RunParam','RunConfig');
    end
end
