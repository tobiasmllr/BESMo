function [ rundata ] = netcdf_import_all( nc_filename, rundata, goThroughSubsurfaceLayers )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% load dimensions and variables
try
    ncid = netcdf.open(nc_filename,'NC_NOWRITE');
catch
    % seperate error, as we might want to just wait for the file to become
    % available... or so... just normal error might prove to complicate
    % things later on if accessing files that are still being written to
    % (which usually should be no problem. Otherwise try to get a copy of a
    % file...
    error('NetCDF error...')
end

%ncdisp(nc_filename)

% get NetCDF dimensions
[dimname0, prints]      = netcdf.inqDim(ncid,0);         % prints (temporal):
[dimname1, nodes_N]     = netcdf.inqDim(ncid,1);        % nodes (spatial)
[dimname2, gsd]         = netcdf.inqDim(ncid,2);            % gsd (number of mean grain size classes, one less as geometric mean)
[dimname3, gsdp1]       = netcdf.inqDim(ncid,3);          % gsd+1 (number of grain size classes)
[dimname4, dim_const]   = netcdf.inqDim(ncid,4);
[dimname5, dim_string]  = netcdf.inqDim(ncid,5);
[dimname6, LstrMat]     = netcdf.inqDim(ncid,6);  % LstrMat (number of subsurface layers


try
    %% get constants
    % metadata
    rundata.nc_current_run         = ncread(nc_filename,'current_run');
    rundata.nc_run_name            = ncread(nc_filename,'run_name');
    rundata.nc_creation_date       = ncread(nc_filename,'creation_date');
    % fixed values (got these from mat-import)
    %     rundata.param_dt            = ncread(nc_filename,'param_dt');
    %     rundata.param_reachlength   = ncread(nc_filename,'reach_length');
    rundata.param_reachbankfull = ncread(nc_filename,'reach_width');
    %     rundata.dx_array            = ncread(nc_filename,'dx_array');
    %     rundata.param_delta_s       = ncread(nc_filename,'param_delta_s');
    %     rundata.zbase               = ncread(nc_filename,'zbase');
    rundata.gsd_initDsi         = ncread(nc_filename,'gsd_initDsi');
catch
    error('NetCDF error...')
end

try
    %% get temporal variables
    % zero dimensional var
    %rundata.store_t       = ncread(nc_filename,'store_t');
    %rundata.store_time    = ncread(nc_filename,'store_time');
    
    % one dimensional var (nodesN)
    rundata.store_etab      = ncread(nc_filename,'store_etab');
    rundata.store_Dsg       = ncread(nc_filename,'store_Dsg');
    rundata.store_Ds90      = ncread(nc_filename,'store_Ds90');
    %rundata.store_SDg       = ncread(nc_filename,'store_SDg');
    rundata.store_qbx       = ncread(nc_filename,'store_qbx');
    %rundata.store_Msj       = ncread(nc_filename,'store_Msj');
    rundata.store_slope     = ncread(nc_filename,'store_slope');
    rundata.store_ustar     = ncread(nc_filename,'store_ustar');
    %rundata.store_waterSurf = ncread(nc_filename,'store_waterSurf');
    rundata.store_TranspPbi = ncread(nc_filename,'store_TranspPbi');
    
    % two dimensional var (nodesN, GSD)
    rundata.store_SurfFsi   = ncread(nc_filename,'store_SurfFsi');
    %rundata.store_SurfPbi   = ncread(nc_filename,'store_SurfPbi');
catch
    error('NetCDF error...')
end

if goThroughSubsurfaceLayers
    %% This data is quite a lot, so only load segments of it...
    rundata.store_Dssg      = NaN(size(rundata.store_Dsg));
    rundata.store_Dss90     = NaN(size(rundata.store_Dsg));
    
    syncInterval  = 50; % number of segmented data loaded=
    store_Dssg_subsurface   = NaN(nodes_N,LstrMat,syncInterval);
    store_Dss90_subsurface  = NaN(nodes_N,LstrMat,syncInterval);
    store_eta_subs          = NaN(nodes_N,LstrMat,syncInterval);
    store_pssi_subs         = NaN(nodes_N,LstrMat,gsd,syncInterval);
    %ActiveLayerN            = NaN(nodes_N,syncInterval);
    
    t_max = length(rundata.store_etab);
    
    for f=1:1:t_max
        if mod(f,syncInterval)==1
            syncIntervalLoad = min(syncInterval,t_max - f);
            try
                %ActiveLayerN(:,1:syncInterval)           = ncread(nc_filename,'store_Msj',       [1 f], [nodes_N syncInterval]);
                %store_waterSurf(:,1:syncInterval)       = ncread(nc_filename,'store_waterSurf', [1 f], [nodes_N syncInterval]);
                %store_etab(:,1:syncInterval)            = ncread(nc_filename,'store_etab',      [1 f], [nodes_N syncInterval]);
                %store_Dsg(:,1:syncInterval)             = ncread(nc_filename,'store_Dsg',       [1 f], [nodes_N syncInterval]);
                
                %highestActive = max(ActiveLayerN(:)) + 1;
                
                %3D
                store_pssi_subs(:,:,:,1:syncIntervalLoad)       = ncread(nc_filename,'store_pssi_subs',         [1 1 1 f], [nodes_N LstrMat gsd syncIntervalLoad]);
                %store_Dssg_subsurface(1:LstrMat,:,1:syncInterval)  = ncread(nc_filename,'store_Dssg_subsurface',  [1 1 f], [nodes_N LstrMat syncIntervalLoad]);
                %store_Dss90_subsurface(1:LstrMat,:,1:syncInterval) = ncread(nc_filename,'store_Dss90_subsurface', [1 1 f], [nodes_N LstrMat syncIntervalLoad]);
                store_eta_subs(:,:,1:syncIntervalLoad)          = ncread(nc_filename,'store_eta_subs',         [1 1 f], [nodes_N LstrMat syncIntervalLoad]);

            catch
                error('NetCDF error...')
            end
            
            for i=1:syncIntervalLoad
                [store_Dssg_subsurface(:,:,i),store_Dss90_subsurface(:,:,i),~,~] = Function_getSubsurfaceDg(store_pssi_subs(:,:,:,i),store_eta_subs(:,:,i),LstrMat,nodes_N,rundata.gsd_initDsi);
            end
            
            sync = 1;
        end
        for n=1:nodes_N
            rundata.store_Dssg(n,f) = store_Dssg_subsurface(n,1,sync);
            rundata.store_Dss90(n,f)= store_Dss90_subsurface(n,1,sync);
        end
        sync = sync + 1;
    end
    
    %% find how many time steps are in the file
    % try
    %     store_time = ncread(nc_filename,'store_time');
    % catch
    %     error('NetCDF error...')
    % end
    % index=sum(~isnan(store_time));
    %
    %
    % nodes_loc = zeros(length(dx_array),1);
    % nodes_loc(2:end) = cumsum(dx_array(1:end-1));
    % nodes_loc2D = transpose(repmat(nodes_loc,[1 LstrMat]));
    
    
end

clear store_Dssg_subsurface store_Dss90_subsurface store_eta_subs store_pssi_subs ActiveLayerN
close all
