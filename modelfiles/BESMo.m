function [ ] = FlumeModel( RunParam, current_run, previous_run, outputStyle, nc_params, run_mode, scratchDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
format long e

if run_mode.pickup && exist(strcat(current_run,'_varDump.mat'),'file')
    disp('Picking up run, no need to load files!')
else
    %% Loading parametrization
    % Loading Physical variables
    %[file.PhysicalData] = load('PhysicalData.txt');
    param_reachlength   = RunParam.ReachLength;         % Length of the reach, in m
    %param_reachslope    = RunParam.Slope0;              % Initial bed Slope, in m/m
    %param_reachbankfull = RunParam.ChBankfullWidth;     % Channel Bankfull width, in m
    %param_feedposition  = RunParam.FeedPosition;
    
    % Loading Physical constants
    %[file.PhysicalConstants] = load('PhysicalConstants.txt');
    const_rhoS      = RunParam.RhoS;            % Sediment density, in kg/m3
    const_lambda    = RunParam.Lambda;          % bed porosity
    const_g         = RunParam.Gravity;         % Acceleration of gravity, in m/s2
    const_rhoW       = RunParam.RhoW;            % Water density, in kg/m3
    const_r         = const_rhoS/const_rhoW - 1;     % Specifc submeged gravity
    
    % Auxiliary parameters
    %[file.AuxParemeters] = load('AuxiliaryParameters.txt');
    param_nk        = RunParam.nk;          % nk: Dimensionless Nikuradse number
    param_alfa_r    = RunParam.alfar;       % alfar: Dimensionless Manning - Strickler coefficient evaluated by means of Parker (1991)
    param_alfa_F    = RunParam.alfaF;       % param_alfa_F: Dimensionless coefficient for exchange fractions between the Active Layer and the Substrate
    param_nLa       = RunParam.nLa;         % nLa: Dimensionless coefficient for the active layer thickness
    param_au        = RunParam.au;          % au: Upwinding coefficient
    
    
    % Loading Time and Space parameters for discretition
    %[file.DiscretData] = load('DiscretData.txt');
    NPrint          = RunParam.NPrint;                  % Number of printouts after the initial one
    param_NtoPrint  = RunParam.NtoPrint;                % Number of steps until a printout is made
    %param_dx        = RunParam.dx;                   	% cell size, in m
    param_dt        = RunParam.dt;                   	% time step, in s
    param_freqQ     = RunParam.freqQ;                 	% Hydrograph frequency
    %param_Thours    = NPrint*param_NtoPrint*param_dt/3600;% hours of calculation time
    stepN           = NPrint*param_NtoPrint;              % total time steps number
    
    % Loading Stratigraphy Data
    %[file.DataStrat] = load('StratigraphyData.txt');
    param_deta_s    = RunParam.deta_s;     % deta_s = Thickness of the layer to store the stratigraphy, in m
    param_alfa_LZu   = RunParam.alfa_LZu;   % alfa_LZu =  Dimensionless parameter that define top bed deposit elevation to consider in the stratigraphy
    param_alfa_LZd   = RunParam.alfa_LZd;   % alfa_LZu =  Dimensionless parameter that define bottom bed deposit elevation to consider in the stratigraphy
    
    
    % reading Incoming  User defined hydrograph
    file.Hyd_Data_Type3 = load('UserDefinedHydrograph.txt');
    param_hydro_Tw        = file.Hyd_Data_Type3(:,1)*3600;
    param_hydro_Qw        = file.Hyd_Data_Type3(:,2);
    
    % net incoming sediment transport rate qin, in m2/s
    file.qin = load('qb_feed.mat');
    % transp_feed_qin_T = file.qin.TArray .* 3600;
    feed_TArray_qsi_const = file.qin.FeedArrayConst;
    feed_TArray_qsi_pulsed = file.qin.FeedArrayPulsed;
    
    clear file
    
    %% Discretization variables
    % number of computing nodes, including the extrems
    %nodes_N = round(param_reachlength/param_dx) + 1;
    % location vector: locations of nodes
    %nodes_loc = 0 : param_dx : param_reachlength;
    % temporal vector: each value corresponds to the staring time of a given time step
    
    dx_array = RunParam.dx_array;
    
    nodes_loc = RunParam.nodes_loc;
    nodes_N = length(nodes_loc);
    %time = (0 : stepN - 1) * param_dt;
    
    %% Loading Sediment transport capacity formula
    %          1: Parker (1990)
    %          2: Wilcock-Crowe (2003)
    %          3: Ashida-Michiue (1972)
    %          4: Meyer-Peter Muller
    
    if RunParam.eqNumber == 1
        % Parker Strain Functions
        %RunParam.Strfunc = load('Function_Parker_Strain.txt');
        parkerStrain.phisg0_ST   = RunParam.Strfunc(:,1);
        parkerStrain.w0_ST       = RunParam.Strfunc(:,2);
        parkerStrain.s0_ST       = RunParam.Strfunc(:,3);
        
        [ parkerStrain.intervect,parkerStrain.omega0inter,parkerStrain.sigma0inter ] ...
            = Calc_ParkerStrain( parkerStrain.phisg0_ST, parkerStrain.w0_ST, parkerStrain.s0_ST );
    end
    
    
    %% Load inital GSDs
    % Surface
    [gsd_initDsiFile,gsd_initSurfFsiFile] ...
        = importGSD('GSD_InitialSurf.txt',2);
    [~,D_index] = unique(gsd_initDsiFile);
    gsd_initDsi = gsd_initDsiFile(sort(D_index));
    gsd_initSurfFsi = gsd_initSurfFsiFile(sort(D_index));
    
    % number of sediment tractions
    gsd_MG = length(gsd_initSurfFsi)-1;
    % convert Fsi from grain size array to spatial array:
    gsd_initSurfFsi = transpose(repmat(gsd_initSurfFsi(:),1,nodes_N));
    
    % calc grain size statistics:
    [gsd_initSurfPfi,gsd_initSurfDg,gsd_initSurfSDg,gsd_initSurfD90,~] ...
        = Calc_GSD_PfiDgSdgD90Dimean_noArray(gsd_initSurfFsi,gsd_initDsi,gsd_MG,nodes_N);
    
    % Subsurface
    [~,gsd_initSubsFsiFile] ...
        = importGSD('GSD_InitialSubSurf.txt',2);
    gsd_initSubsFsi = gsd_initSubsFsiFile(sort(D_index));
    % convert Fsi from grain size array to spatial array:
    gsd_initSubsFsi = transpose(repmat(gsd_initSubsFsi(:),1,nodes_N));
    % calc grain size statistics:
    [gsd_initSubsPfi,~,~,~,gsd_initDimean] ...
        = Calc_GSD_PfiDgSdgD90Dimean_noArray(gsd_initSubsFsi,gsd_initDsi,gsd_MG,nodes_N);
    
    % Feed
    [gsd_feedDsiFile,gsd_feedFsiFile] ...
        = importGSD('GSD_pbfeed.txt',2);
    gsd_feedDsi = gsd_feedDsiFile(sort(D_index));
    gsd_feedFsi = gsd_feedFsiFile(sort(D_index));
    
    % convert Fsi from grain size array to spatial array:
    gsd_feedFsi = transpose(repmat(gsd_feedFsi(:),1,nodes_N));
    % calc grain size statistics:
    [gsd_feedPfi_Arr,~,~,~,~] ...
        = Calc_GSD_PfiDgSdgD90Dimean_noArray(gsd_feedFsi,gsd_feedDsi,gsd_MG,nodes_N);
    gsd_feedPfi = gsd_feedPfi_Arr(1,:);
    
    
    % Location of the points where saving the frequencies
    %gsd_xGSD = load('xGSD_evolution.txt');
    %gsd_MGSD = length(gsd_xGSD);
    
    clear file;
    
    % =========================================================================
    % DONE INIT
    % =========================================================================
    
    %% pre-alloating variables. Make the code faster
    NodeSubs.ActiveLayerZ_old = zeros(nodes_N,1);    %Active Layer Thickness
    NodeSubs.ActiveLayerZ = zeros(nodes_N,1);        %Active Layer Thickness
    
    ActiveLayerN = zeros(nodes_N,1);        %Active Layer Number
    % Surface    GSD array Variables (spatial arrays)
    NodeGSD.SurfDg = zeros(nodes_N,1);
    NodeGSD.SurfSDg = zeros(nodes_N,1);
    NodeGSD.SurfD90 = zeros(nodes_N,1);
    NodeGSD.SurfD90_new = zeros(nodes_N,1);
    NodeSedTrans.qb = zeros(nodes_N,1);      % Transport rate:
    
    %% Load initial values
    % Initial Surface GS Properties, all the same for all location nodes
    NodeGSD.SurfDg(:) = gsd_initSurfDg;
    NodeGSD.SurfSDg(:) = gsd_initSurfSDg;
    NodeGSD.SurfD90(:) = gsd_initSurfD90;
    NodeGSD.SurfPfi = gsd_initSurfPfi;
    NodeGSD.SurfFsi = gsd_initSurfFsi;
    
    %% Linear interpolation of the incoming water discharge
    % Note: the final time of the hydrograph must be, at least, as large as the
    % time at which computation finishes.
    Lw = length(param_hydro_Qw);
    mi = zeros(Lw-1,1);
    Qw0 = zeros(Lw-1,1);
    Tw0 = zeros(Lw-1,1);
    
    for n = 1:Lw-1
        mi(n,1) = (param_hydro_Qw(n+1,1) - param_hydro_Qw(n,1))/((param_hydro_Tw(n+1,1) - param_hydro_Tw(n,1)));
        Qw0(n,1) = param_hydro_Qw(n,1);
        Tw0(n,1) = param_hydro_Tw(n,1);
    end
    
    
    %% Initial Conditions (bed and Surface GSD)
    
    % Initial bed elevation, in m
    NodeGeom.Elev = RunParam.elev_array;
    NodeGeom.Slope = RunParam.slope_array;
    
    
    % Initial stratigrahy
    zbase = NodeGeom.Elev(nodes_N,1) - param_alfa_LZd*NodeGeom.Elev(1,1);
    ztop = (1+param_alfa_LZu)*NodeGeom.Elev(1,1);
    
    % number of storing points within the deposit at each computational node.
    % Active layer is excluded
    for n = 1:nodes_N
        La_str = param_nLa*NodeGSD.SurfD90(n,1)/1000;
        Msi = floor((NodeGeom.Elev(n,1)-La_str-zbase)/param_deta_s) + 2;
        ActiveLayerN(n,1) = Msi;
    end
    
    % Number of storing points. Two extra-points are added to account for the
    % bed elevation and the active layer-substrate interphase regardless they
    % fall in the deta_s spacing
    if rem((ztop - zbase),param_deta_s) == 0
        LstrMat = round((ztop - zbase)/param_deta_s) + 1 + 2;
    else
        LstrMat = floor((ztop - zbase)/param_deta_s + 2) + 2;
    end
    %% Prepare store variables
    % zero dimensional var
    store.t             = NaN(nc_params.syncInterval, 1);
    store.time          = NaN(nc_params.syncInterval, 1);
    
    % one dimensional var (nodesN)
    store.BedElev       = NaN(nc_params.syncInterval, nodes_N);
    store.SurfDg        = NaN(nc_params.syncInterval, nodes_N);
    store.SurfD90       = NaN(nc_params.syncInterval, nodes_N);
    store.SurfSDg       = NaN(nc_params.syncInterval, nodes_N);
    store.SurfQbx       = NaN(nc_params.syncInterval, nodes_N);
    store.ActiveLayerN  = NaN(nc_params.syncInterval, nodes_N);
    store.slope         = NaN(nc_params.syncInterval, nodes_N);
    store.ustar         = NaN(nc_params.syncInterval, nodes_N);
    store.waterd        = NaN(nc_params.syncInterval, nodes_N);
    
    % two dimensional var (nodesN, GSD)
    store.TranspQbi     = NaN(nc_params.syncInterval, nodes_N, gsd_MG);
    store.SurfFsi       = NaN(nc_params.syncInterval, nodes_N, gsd_MG+1);
    store.SurfPfi       = NaN(nc_params.syncInterval, nodes_N, gsd_MG);
    
    % two dimensional var (nodesN, LStrMat)
    store.SubsEta       = NaN(nc_params.syncInterval, nodes_N, LstrMat);
    
    % three dimensional var (nodesN, LstrMat, GSD)
    store.SubsPssi      = NaN(nc_params.syncInterval, nodes_N, LstrMat, gsd_MG);
end

%% PREPARE STORAGE TO FILES
%% Create storage .mat file
%storemat_filename = strcat(current_run,'_store.mat');
%save(storemat_filename,'current_run','-v7.3');
%storemat = matfile(storemat_filename,'Writable',true);

%fprintf(fdev,strcat('Store Variables as: ',storemat_filename,'\n'));
%disp(strcat('Store Variables as: ',storemat_filename))

%% STORE DATA IN NETCDF
if run_mode.pickup && exist(strcat(current_run,'_varDump.mat'),'file')
    % load previous nc
    % populate initial variables
    % set current time-step
    
    % load snapshot of variables in workspace:
    load(strcat(current_run,'_varDump.mat'))
    
    assert(synci==1,'synci was not equal 1 when the snapshot was saved! This implies the snapshot was taken in the wrong state.')
    warning('Picking up previous run. Experimental feature!')
    
    % advance to step after save:
    index = index+1;
    t_before = t_step;
    fprintf('LOAD SIMULATION at t = %i, index = %i... ',t_step,index);
    
    % load netcdf:
    nc_filename = strcat(current_run,'_store.nc');
    if exist(nc_filename,'file')
        try
            ncid = netcdf.open(nc_filename,'WRITE');
        catch err
            fprintf('Could not reopen netcdf file....')
        end
    end
    
    %% next iteration
    store.qb_mean = 0 .* NodeSedTrans.qb;
    store.pbi_mean = 0 .* NodeSedTrans.qbi;
    t_step = t_step + 1;
    tic
elseif run_mode.pickup && (exist(strcat(scratchDir,'/',previous_run,'_varDump.mat'),'file') || exist(strcat(previous_run,'_varDump.mat'),'file'))
    %% pick up a run that had the same parameters, but was shorter than the new one
    % load data (Overwrites all data with the old stuff)
    % use regexp to exclude the variables that start with the listed
    % names, as these should not be overwritten with the old parameters
    if exist(strcat(scratchDir,'/',previous_run,'_varDump.mat'),'file')
        load(strcat(scratchDir,'/',previous_run,'_varDump.mat'),'-regexp', '^(?!NPrint|stepN|current_run|feed_TArray_qsi_const|feed_TArray_qsi_pulsed)...')
    else
        load(strcat(previous_run,'_varDump.mat'),'-regexp', '^(?!NPrint|stepN|current_run|feed_TArray_qsi_const|feed_TArray_qsi_pulsed)...')
    end
    
    assert(synci==1,'synci was not equal 1 when the snapshot was saved! This implies the snapshot was taken in the wrong state.')
    warning('Picking up previous run. Experimental feature!')
    
    t_step = time;
    % advance to step after save:
    index = index+1;
    t_before = t_step;
    fprintf('LOAD SIMULATION at t = %i, index = %i... ',t_step,index);
    
    store.current_run               = current_run;
    store.run_name                  = current_run;
    
    % we have to create a new netcdf file
    if nc_params.saveNC_subs || nc_params.saveNC
        [ ncid, nc_obj,nc_params ] = netcdf_create( current_run, NPrint,nodes_N,gsd_MG,LstrMat, nc_params );
    end
    
    % load old netcdf:
    nc_filename_previous = strcat(current_run,'_store.nc');
    if exist(nc_filename_previous,'file')
        try
            ncid_previous = netcdf.open(nc_filename_previous,'WRITE');
        catch err
            fprintf('Could not reopen previous netcdf file....')
        end
    end
    
    
    
    %% next iteration
    t_step = t_step + 1;
    tic
else
    %% START A NEW RUN
    if nc_params.saveNC_subs || nc_params.saveNC
        [ ncid, nc_obj,nc_params ] = netcdf_create( current_run, NPrint,nodes_N,gsd_MG,LstrMat, nc_params );
    end
    %% Initial stratigraphy
    [~,NodeSubs.Elev,NodeSubs.Pfi,pssi] = InitialStratigraphy(LstrMat,nodes_N,gsd_MG,...
        nodes_loc,param_nLa,NodeGSD.SurfD90,ActiveLayerN,NodeGeom.Elev,...
        param_deta_s,zbase,gsd_initSubsPfi,NodeGSD.SurfPfi);
    
    
    RunParam.nLa = param_nLa;
    RunParam.detas = param_deta_s;
    % switch old geometry to new one:
    for n=1:nodes_N
        shdist = ActiveLayerN(n);
        NodeSubs.Elev(n,:) = circshift(flip(NodeSubs.Elev(n,:),2),shdist+1,2);
        NodeSubs.Pfi(n,:,:) = circshift(flip(NodeSubs.Pfi(n,:,:),2),shdist+1,2);
        
        % fill all bottom area:
        fidx = find(isnan(NodeSubs.Elev(n,:)),1,'first');
        lowerval = NodeSubs.Elev(n,fidx(1)-1);
        NodeSubs.Elev(n,fidx:end) = [-RunParam.detas:-RunParam.detas:-RunParam.detas*(LstrMat-fidx+1)] + lowerval;
    end
    
    % Calculation of the active layer thickness, La (in m)
    NodeSubs.ActiveLayerZ_old(:) = 0;
    NodeSubs.ActiveLayerZ(:) = param_nLa*NodeGSD.SurfD90(:)/1000;
        
    [NodeGSD.SurfDg(:),NodeGSD.SurfSDg(:),NodeGSD.SurfD90_new(:),NodeGSD.SurfD50,~] = ...
        GSDComputing_vectorized_spatial(NodeGSD.SurfFsi(:,:),gsd_initDsi,gsd_MG,nodes_N);
    
    NodeGeom.WidthFixed = RunParam.ChBankfullWidth_array;
    NodeGeom.WidthFlow = RunParam.ChBankfullWidth_array;
    
    time = 0;
    for n = 1:Lw-1
        if (time >= param_hydro_Tw(n,1)) && (time < param_hydro_Tw(n+1,1))
            qwi_i = (Qw0(n,1) + mi(n,1)*(time - Tw0(n,1)))/RunParam.ChBankfullWidth_array(n);
            qwi_s = (Qw0(n,1) + mi(n,1)*(time + param_dt - Tw0(n,1)))/RunParam.ChBankfullWidth_array(n);
            hydro_TArray_qwi = 0.5*(qwi_i + qwi_s);
        end
    end
    
    [ustar,waterDepth,~] = Function_BackwaterCalc(hydro_TArray_qwi,...
        NodeGeom.Slope,NodeGeom.Elev, NodeGSD.SurfD90_new,...
        dx_array,nodes_loc,param_nk,param_alfa_r,...
        const_rhoW,const_g,RunParam.backwaterparam);
    
    % Calculate initial bedload frequencies:
    switch RunParam.eqNumber
        case 1
            % Parker's sediment transport function
            %[~,temp_pbqs] = Function_Parker(gsd_initDsi,gsd_initSurfPfi,NodeGSD.SurfDg(1),NodeGSD.SurfSDg(1),const_r,const_g,ustar,parkerStrain.intervect,parkerStrain.omega0inter,parkerStrain.sigma0inter);
            [~,NodeSedTrans.qbi] = Function_Parker_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi,NodeGSD.SurfDg,NodeGSD.SurfSDg,const_r,const_g,ustar,parkerStrain.intervect,parkerStrain.omega0inter,parkerStrain.sigma0inter,gsd_MG,nodes_N,0);
        case 2
            % Wilcock and Crowe's sediment transport function
            [~,NodeSedTrans.qbi] = Function_WilcockCrowe_vectorized_spatial(gsd_initDsi,gsd_initDimean,NodeGSD.SurfPfi,NodeGSD.SurfDg,const_r,const_g,ustar,gsd_MG,nodes_N,0,RunParam.tau_crit_factor);
        case 3
            % Ashida and Michiue's sediment transport function
            [~,NodeSedTrans.qbi] = Function_AshidaMichiue1972_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi,NodeGSD.SurfDg,const_r,const_g,ustar,gsd_MG,nodes_N,0);
        case 4
            % Meyer-Peter Muller
            [~,NodeSedTrans.qbi] = Function_MPM_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi,NodeGSD.SurfDg,const_r,const_g,ustar,gsd_MG,nodes_N,0,RunParam.tau_crit_factor);
        otherwise
            error('eqNumber to choose sediment transport function out of bound')
    end
    
    %% Saving initial data
    store.qb_mean = NodeSedTrans.qb;
    store.pbi_mean = NodeSedTrans.qbi;
    index = 0;
    synci = 1;
    storei = 1;
    t_step = 0;
    
    if nc_params.saveNC_subs || nc_params.saveNC
        % zero dimensional var
        store.t(synci)                  = t_step;
        store.time(synci)               = t_step*param_dt/3600;
        
        % one dimensional var (nodesN)
        store.BedElev(synci,:)          = NodeGeom.Elev;
        store.SurfDg(synci,:)           = NodeGSD.SurfDg;
        store.SurfD90(synci,:)          = NodeGSD.SurfD90;
        store.SurfSDg(synci,:)          = NodeGSD.SurfSDg;
        store.SurfQbx(synci,:)          = NodeSedTrans.qb;
        store.ActiveLayerN(synci,:)     = ActiveLayerN;
        store.slope(synci,:)            = NodeGeom.Slope;
        store.ustar(synci,:)            = ustar;
        store.waterd(synci,:)           = waterDepth;
        
        % two dimensional var (nodesN, GSD)
        store.TranspQbi(synci,:,:)      = NodeSedTrans.qbi;
        store.SurfFsi(synci,:,:)        = NodeGSD.SurfFsi;
        store.SurfPfi(synci,:,:)        = NodeGSD.SurfPfi;
        %store.SubsDg(synci,:,:)         = gsd_array_SubsDg_save;
        %store.SubsD90(synci,:,:)        = gsd_array_SubsD90_save;
        store.SubsEta(synci,:,:)        = NodeSubs.Elev;
        
        % two dimensional var (nodesN, LStrMat)
        %store.SubsEta       = NaN(nc_params.syncInterval, LstrMat, nodes_N);
        
        % three dimensional var (nodesN, LstrMat, GSD)
        store.SubsPssi(synci,:,:,:)      = NodeSubs.Pfi;
                
        % OTHER NON TEMPORAL, saved only once
        store.current_run               = current_run;
        store.run_name                  = current_run;
        store.dx_array                  = dx_array;
        store.param_dt                  = param_dt;
        store.param_reachlength         = param_reachlength;
        store.reachwidth                = RunParam.ChBankfullWidth_array;
        store.gsd_initDsi               = gsd_initDsi;
    end
    index = index + 1;
    
    if outputStyle.showStoreSize
        ncfile = dir(nc_params.filename);
        ncfile.sizeMB = ncfile.bytes/1024/1024;
        message = sprintf('Netcdf file size in MB: %.2f, with deflate level: %i',ncfile.sizeMB,nc_params.deflate_level);
        disp(message)
    end
    
    %% temporal loop
    
    
    doingSubtimeStep = false;
    hadNegativeSlopes = false;
    %subtimestepsArray = ones(stepN,1);
    %overThresholdFactorArrayMax = NaN(stepN,1);
    %overThresholdFactorArrayAfterRecalc = NaN(stepN,1);
    
    %CourantNumberArray = NaN(stepN,1);
    
    message = sprintf('Progress in Run: %02d %%\n',0);
    disp(message)
    
    tic
    t_step = 1;
    t_before = 0; % just for checking repeating time steps for storing values
    subtimesteps = 1;
    output_pfeed_start = false;
end

while t_step < (stepN + 1)
    %% OUTPUT FOR LOGFILE
    if (mod(t_step,round(stepN/100)) == 0)
        percentage = round((t_step/stepN)*100);
        secPerStep = toc;
        minExpected = round((100 - percentage) * secPerStep / 60);
        hourExpected = floor(minExpected / 60);
        minAfterHourExpected = minExpected - (hourExpected * 60);
        message = sprintf('Progress in Run: %02d %%, time left (hour:min): %02d:%02d, 1 percent took: %d sec\n',percentage,hourExpected,minAfterHourExpected, round(secPerStep));
        disp(message)
        tic
    end
    
    RunParam.subtimestep_max = RunParam.dt * 25;
    RunParam.CourantThresh = 0.1;
    RunParam.OnlyCourantThresh = true;
    
    if run_mode.doSubtimesteps && t_step>1
        [t_step, tt, t_day,t_date, index, subtimesteps, RunParam.dt_actual , ...
            NodeSedTrans,NodeSubs,NodeGeom,NodeGSD,doingSubtimeStep,...
            subtimestepsArray,overThresholdFactorArrayMax, ...
            overThresholdFactorArrayAfterRecalc, ...
            CourantNumber,increaseSubtimestep] = Function_dynamicTimestepThreshold_LosPadres(...
            nodes_N, t_step, RunParam, outputStyle, ...
            NodeSedTrans,NodeSubs,NodeGeom,NodeGSD,timeOutput,SCNode, ...
            subtimestepsArray,overThresholdFactorArrayMax,...
            overThresholdFactorArrayAfterRecalc,...
            subtimesteps,doingSubtimeStep,...
            NodeSedTrans_previous,NodeSubs_previous,NodeGeom_previous,NodeGSD_previous,...
            t_day,t_date,index,t_day_previous,t_date_previous,index_previous,increaseSubtimestep);
        error('subtimesteps removed right now')
    else
        tt=1;
        param_dt_actual=param_dt;
    end
    
    CourantNumber = max(NodeSedTrans.qb) * param_dt / min(dx_array);
    
    if CourantNumber > 0.5
        warning('Courant Number > 0.5. Change timestep or spatial resolution!')
    end
    
    % Discharge entering the flume:
    time = t_step;
    for n = 1:Lw-1
        if (time >= param_hydro_Tw(n,1)) && (time < param_hydro_Tw(n+1,1))
            qwi_i = (Qw0(n,1) + mi(n,1)*(time - Tw0(n,1)))/RunParam.ChBankfullWidth_array(n);
            qwi_s = (Qw0(n,1) + mi(n,1)*(time + param_dt - Tw0(n,1)))/RunParam.ChBankfullWidth_array(n);
            hydro_TArray_qwi = 0.5*(qwi_i + qwi_s);
        end
    end
    qw = hydro_TArray_qwi;
    
    % Sediment entering the flume:
    feed_array_qin = zeros(nodes_N,1);
    feed_array_qin(RunParam.constfeednode)  = feed_TArray_qsi_const(t_step);
    
    % spread pulsed feed over N nodes:
    feednodes=1; % can't be even number
    if feednodes > 1
        feed_array_qin(RunParam.pulsedfeednode-floor(feednodes/2):RunParam.pulsedfeednode+floor(feednodes/2)) = feed_TArray_qsi_pulsed(t_step)/feednodes + feed_array_qin(RunParam.pulsedfeednode-floor(feednodes/2):RunParam.pulsedfeednode+floor(feednodes/2));
        %feed_array_qin(RunParam.constfeednode) = feed_TArray_qsi_pulsed(t);
    else
        % If the pulsed feed is entering only one node, then add the pulsed
        % feed to what already is fed here (the pulsed feed node can be the
        % same as the constant feed node.
        feed_array_qin(RunParam.pulsedfeednode) = feed_array_qin(RunParam.pulsedfeednode) + feed_TArray_qsi_pulsed(t_step);
    end
    
    if (feed_TArray_qsi_pulsed(t_step) > 0) && (output_pfeed_start == false)
        output_pfeed_start  = true;
        fprintf('%s ## index: %i, Feed START at step=%i, hour=%.2f\n',current_run,index,t_step, t_step*param_dt/3600)
    end
    
    if (feed_TArray_qsi_pulsed(t_step) == 0) && (output_pfeed_start == true)
        output_pfeed_start  = false;
        fprintf('%s ## index: %i, Feed END at step=%i, hour=%.2f\n',current_run,index,t_step, t_step*param_dt/3600)
    end
    
    % store the current sediment transport (in case we have to repeat
    % timestep)
    %sedTrans_array_qbx_previous = NodeSedTrans.qb;
    
    while tt < (subtimesteps + 1)
        timeOutput = (t_step - 1) + (tt/subtimesteps); % Only used for display
        
        %% USUAL LOOP:
        
        %% Calculate new Bed Elevation and difference in transport rate:
        %         [NodeGeom.Elev_new, sedTrans_array_dqbx, NodeGeom.dElev_feed] = Function_GetNewElevation4( ...
        %             NodeGeom.Elev, NodeSedTrans.qb, feed_array_qin, nodes_N, dx_array, ...
        %             param_dt_actual, param_au, const_lambda, param_freqQ, timeOutput,RunParam.ChBankfullWidth_array);
        %
        
        NodeSedTrans.qb = NodeSedTrans.qb;
        NodeSedTrans.qb_feed= feed_array_qin;
        %NodeGeom.Elev   = NodeGeom.Elev;
        NodeGeom.dx     = dx_array;
        RunParam.Upwind_au  = param_au;
        RunParam.Upwind_order = 'first';
        RunParam.Lambda     = const_lambda;
        RunParam.dt_actual  = param_dt_actual;
        RunParam.freqQ      = param_freqQ;
        
        [NodeGeom, NodeSedTrans] = Function_GetNewElevation( ...
            NodeGeom, NodeSedTrans, nodes_N, RunParam, timeOutput);
        
        %sedTrans_array_dqbx = NodeSedTrans.dqb;
                
        NodeGSD.FeedPfi = transpose(repmat(gsd_feedPfi(:),1,nodes_N));
        
        [ NodeGSD,NodeSedTrans ] = Function_GetNewParticleSizeDiff( ...
                    NodeGeom, NodeGSD, NodeSedTrans, NodeSubs,RunParam,t_step );
                
        %% Recalc statistics vectorized
        [NodeGSD.SurfDg(:),NodeGSD.SurfSDg(:),NodeGSD.SurfD90_new(:),NodeGSD.SurfD50,NodeGSD.SurfPfi(:,:)] = ...
            GSDComputing_vectorized_spatial(NodeGSD.SurfFsi(:,:),gsd_initDsi,gsd_MG,nodes_N);
        
        
        assert(all(NodeGSD.SurfDg<gsd_initDsi(end)),'error with GS');
        NodeGSD.SurfPfi_new = NodeGSD.SurfPfi;
        
        %Error if problem
        if isnan(sum(NodeGSD.SurfD90_new(:)))
            errorposition   = 'Error:  j = %i, i = %i \n';
            errorcause      = 'isnan(NodeGSD.SurfD90_new(:,1))';
            error(strcat(errorposition,errorcause), timeOutput, n)
        end
        
        %% Calculation of the active layer thickness, La (in m)
        NodeSubs.ActiveLayerZ_old(:) = NodeSubs.ActiveLayerZ(:);
        NodeSubs.ActiveLayerZ(:) = param_nLa*NodeGSD.SurfD90_new(:)/1000;
        
        Reach.NodeToReach = [1:length(NodeGeom.Elev)]; % all nodes
        RunParam.LagoonStage_fixed = false;
        
        [ustar,waterDepth] = Function_BackwaterCalc3(qw,1,1,...
            NodeGeom, NodeGSD, RunParam, Reach.NodeToReach);
        %         ustar_new = ustar_array;
        %         waterDepth_new = waterDepth;
        %         newtime = newtime + toc;
        %         %
        %         %if rem(t,param_NtoPrint)==0
        %         figure;
        %         plot(NodeGeom.Elev_new,'-k');
        %         hold on;
        %         plot(NodeGeom.Elev_new+waterDepth_old,'-ob');
        %         hold on;
        %         plot(NodeGeom.Elev_new+waterDepth_new,'-*b');
        %         hold on;
        %         plot(ustar_old,'-or');
        %         hold on;
        %         plot(ustar_new,'-*r');
        %         hold off;
        
        %end
        
        %% Calculate Bedload transport
        switch RunParam.eqNumber
            case 1
                % Parker's sediment transport function
                [NodeSedTrans.qb(:),NodeSedTrans.qbi(:,:)] = Function_Parker_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi_new(:,:),NodeGSD.SurfDg(:),NodeGSD.SurfSDg(:),const_r,const_g,ustar,parkerStrain.intervect,parkerStrain.omega0inter,parkerStrain.sigma0inter,gsd_MG,nodes_N,timeOutput);
            case 2
                % Wilcock and Crowe's sediment transport function
                [NodeSedTrans.qb(:),NodeSedTrans.qbi(:,:)] = Function_WilcockCrowe_vectorized_spatial(gsd_initDsi,gsd_initDimean,NodeGSD.SurfPfi_new(:,:),NodeGSD.SurfDg(:),const_r,const_g,ustar(:),gsd_MG,nodes_N,timeOutput,RunParam.tau_crit_factor);
                %                                                 Function_WilcockCrowe_vectorized_spatial(Di,         Dimean,        piN,                 DgN,                    R,            g,ustar,         gsd_MG,nodes_N)
            case 3
                % Ashida and Michiue's sediment transport function
                [NodeSedTrans.qb(:),NodeSedTrans.qbi(:,:)] = Function_AshidaMichiue1972_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi_new(:,:),NodeGSD.SurfDg(:),const_r,const_g,ustar(:),gsd_MG,nodes_N,timeOutput);
                %                                                 Function_AshidaMichiue1972_vectorized_spatial(Dimean,        piN,                  DgN,                   R,      g,      ustar,         gsd_MG,nodes_N)
            case 4
                % Meyer-Peter Muller
                [NodeSedTrans.qb(:),NodeSedTrans.qbi(:,:)] = Function_MPM_vectorized_spatial(gsd_initDimean,NodeGSD.SurfPfi_new(:,:),NodeGSD.SurfDg(:),const_r,const_g,ustar(:),gsd_MG,nodes_N,timeOutput,RunParam.tau_crit_factor);
            otherwise
                error('eqNumber to choose sediment transport function out of bound')
        end
        
        %% Obtaining data for the substrate
        
        [~,NodeSubs] = SubstrateGSD_vectorized(...
            RunParam,NodeGSD,NodeGeom,NodeSubs,...
            NodeSedTrans.fIi,LstrMat,gsd_MG);
        
        %% Restablishing values of eta, psi, pssi, GSD parameters for the loop
        % and the values of the stratigrahy
        
        NodeGeom.Elev = NodeGeom.Elev_new;
        NodeGSD.SurfPfi = NodeGSD.SurfPfi_new;
        %NodeGSD.SurfDg = NodeGSD.SurfDg_new;
        %NodeGSD.SurfSDg = NodeGSD.SurfSDg;
        NodeGSD.SurfD90 = NodeGSD.SurfD90_new;
        %ActiveLayerN = ActiveLayerN_oldVersion;
        
        % END OF SUBTIMESTEP
        tt = tt + 1;
    end
    
    
    %% Saving results
    % keep track of mean transport rates:
    store.qb_mean = store.qb_mean + NodeSedTrans.qb;
    store.pbi_mean = store.pbi_mean + NodeSedTrans.qbi;
    
    if rem(t_step,param_NtoPrint)==0
        if outputStyle.dispSubtime
            %fprintf('%s ## index: %i\n',current_run,index);
            fprintf('%s ## index: %i, simhour: %.2f, progperc = %2.2f\n',current_run,index,t_step*param_dt/3600,(t_step/stepN)*100);
        end
        %figure;plot(NodeGeom.Elev_new,'k-');hold on;plot(NodeGeom.Elev_new + waterDepth,'b-')
        
        %% Stuff we don't have to do every subtimestep (prepare saving things)
        
        if nc_params.saveNC_subs || nc_params.saveNC
            if t_step == t_before
                % if we did repeat this timestep
                % WE ALREADY SAVED THE WRONG VALUE FROM LAST INDEX
                if synci == 1
                    % we synced the step before...
                    synci = nc_params.syncInterval;
                else
                    synci = synci - 1;
                end
            end
            
            % zero dimensional var
            store.t(synci)                  = t_step;
            store.time(synci)               = t_step*param_dt/3600;
            
            % one dimensional var (nodesN)
            store.BedElev(synci,:)          = NodeGeom.Elev;
            store.SurfDg(synci,:)           = NodeGSD.SurfDg;
            store.SurfD90(synci,:)          = NodeGSD.SurfD90;
            store.SurfSDg(synci,:)          = NodeGSD.SurfSDg;
            store.SurfQbx(synci,:)          = store.qb_mean ./ param_NtoPrint;
            %             store.ActiveLayerN(synci,:)     = ActiveLayerN;
            store.slope(synci,:)            = NodeGeom.Slope;
            store.ustar(synci,:)            = ustar;
            store.waterd(synci,:)           = waterDepth;
            
            % two dimensional var (nodesN, GSD)
            store.TranspQbi(synci,:,:)      = store.pbi_mean ./ param_NtoPrint;
            store.SurfFsi(synci,:,:)        = NodeGSD.SurfFsi;
            store.SurfPfi(synci,:,:)        = NodeGSD.SurfPfi;
            %store.SubsDg(synci,:,:)         = gsd_array_SubsDg_save;
            %store.SubsD90(synci,:,:)        = gsd_array_SubsD90_save;
            
            % two dimensional var (nodesN, LStrMat)
            store.SubsEta(synci,:,:)        = NodeSubs.Elev;
            
            % three dimensional var (nodesN, LstrMat, GSD)
            store.SubsPssi(synci,:,:,:)      = NodeSubs.Pfi;
            
            synci = synci + 1;
        end
        
        
        if (rem(index,nc_params.syncInterval) == 0) || (t_step == stepN)
            % sync netcdf to disk every n steps
            if t_step == t_before
                % if we did repeat this timestep
                % WE ALREADY SAVED THE WRONG VALUE FROM LAST INDEX
                storei = storei - 1;
            end
            
            if ((t_step == stepN) && (rem(index,nc_params.syncInterval) ~= 0))
                % if at end of sim and not full store to netcdf:
                restStoreN = rem(index,nc_params.syncInterval);
                
                % zero dimensional var
                restStore.t(1:restStoreN)                  = store.t(1:restStoreN);
                restStore.time(1:restStoreN)               = store.time(1:restStoreN);
                
                % one dimensional var (nodesN)
                restStore.BedElev(1:restStoreN,:)          = store.BedElev(1:restStoreN,:);
                restStore.SurfDg(1:restStoreN,:)           = store.SurfDg(1:restStoreN,:);
                restStore.SurfD90(1:restStoreN,:)          =  store.SurfD90(1:restStoreN,:);
                restStore.SurfSDg(1:restStoreN,:)          = store.SurfSDg(1:restStoreN,:);
                restStore.SurfQbx(1:restStoreN,:)          = store.SurfQbx(1:restStoreN,:);
                restStore.ActiveLayerN(1:restStoreN,:)     = store.ActiveLayerN(1:restStoreN,:);
                restStore.slope(1:restStoreN,:)            = store.slope(1:restStoreN,:);
                restStore.ustar(1:restStoreN,:)            = store.ustar(1:restStoreN,:);
                restStore.waterd(1:restStoreN,:)           = store.waterd(1:restStoreN,:);
                
                % two dimensional var (nodesN, GSD)
                restStore.TranspQbi(1:restStoreN,:,:)      = store.TranspQbi(1:restStoreN,:,:);
                restStore.SurfFsi(1:restStoreN,:,:)        = store.SurfFsi(1:restStoreN,:,:);
                restStore.SurfPfi(1:restStoreN,:,:)        = store.SurfPfi(1:restStoreN,:,:);
                %restStore.SubsDg(1:restStoreN,:,:)         = gsd_array_SubsDg_save;
                %restStore.SubsD90(1:restStoreN,:,:)        = gsd_array_SubsD90_save;
                
                % two dimensional var (nodesN, LStrMat)
                restStore.SubsEta(1:restStoreN,:,:)        = store.SubsEta(1:restStoreN,:,:);
                
                % three dimensional var (nodesN, LstrMat, GSD)
                restStore.SubsPssi(1:restStoreN,:,:,:)      = store.SubsPssi(1:restStoreN,:,:,:);
                
                netcdf_store( ncid, nc_obj, restStore, nc_params, storei, restStoreN )
            else
                % normal save
                netcdf_store( ncid, nc_obj, store, nc_params, storei, nc_params.syncInterval )
            end
            
            netcdf.sync(ncid)
            
            if outputStyle.showStoreSize
                ncfile = dir(nc_params.filename);
                ncfile.sizeMB = ncfile.bytes/1024/1024;
                fprintf('%s ## index: %i, simhour: %.2f, progress = %2.2f\n',current_run,index,t_step*param_dt/3600,round((t_step/stepN)*100));
                fprintf('%s ## Saved NetCDF file. Size in MB: %.2f, deflate level: %i\n',current_run,ncfile.sizeMB,nc_params.deflate_level);
            end
            
            
            synci = 1;
            storei = storei + 1;
            
            if rem(storei,2) == 0
                % Copy the NETCDF file after 50th save to scratch dir, to be able to restart in case...:
                copyfile(ncfile.name,scratchDir);
                
                % take a snapshot of variables in workspace:
                save(strcat(current_run,'_varDump.mat'))
                copyfile(strcat(current_run,'_varDump.mat'),scratchDir);
            end
        end
        
        
        store.qb_mean = 0 .* NodeSedTrans.qb;
        store.pbi_mean = 0 .* NodeSedTrans.qbi;
        
        index = index+1;
        t_before = t_step;
    end
    if run_mode.debug && (run_mode.debugSteps == index)
        fprintf('EXIT SIMULATION at t = %i, index = %i... DEBUG MODE ON',t_step,index);
        return
    end
    %% next iteration
    
    % END OF TIMESTEP
    t_step = t_step + 1;
    
    
    if t_step == (stepN + 1)
        save(strcat(current_run,'.mat'))
        message = sprintf('Returning from FlumeModel loop: Variable snapshot saved as: %s.mat\n',current_run);
        disp(message)
        return
    end
end

%% Save Results
save(strcat(current_run,'.mat'))
message = sprintf('Variable snapshot saved as: %s.mat\n',current_run);
disp(message)
end