function [ ] = netcdf_store( ncid, nc_obj, store, nc_params, storei, storelength )
%NETCDF_STORE Summary of this function goes here
%   Detailed explanation goes here

starti  = (storei - 1) * nc_params.syncInterval;

nodes_N     = size(store.TranspQbi,2);
gsd_MG      = size(store.TranspQbi,3);
LstrMat     = size(store.SubsEta,3);

% zero dimensional var
netcdf.putVar(ncid,nc_obj.var_t           ,[starti], [storelength] ,store.t)
netcdf.putVar(ncid,nc_obj.var_time        ,[starti], [storelength] ,store.time)

% one dimensional var (nodesN)
netcdf.putVar(ncid,nc_obj.var_etab        ,[0 starti],[nodes_N storelength],permute(store.BedElev,[2 1]))
netcdf.putVar(ncid,nc_obj.var_Dsg         ,[0 starti],[nodes_N storelength],permute(store.SurfDg,[2 1]))
netcdf.putVar(ncid,nc_obj.var_Ds90        ,[0 starti],[nodes_N storelength],permute(store.SurfD90,[2 1]))
netcdf.putVar(ncid,nc_obj.var_qbx         ,[0 starti],[nodes_N storelength],permute(store.SurfQbx,[2 1]))
%netcdf.putVar(ncid,nc_obj.var_Msj         ,[0 starti],[nodes_N storelength],permute(store.ActiveLayerN,[2 1]))
netcdf.putVar(ncid,nc_obj.var_slope       ,[0 starti],[nodes_N storelength],permute(store.slope,[2 1]))
netcdf.putVar(ncid,nc_obj.var_ustar       ,[0 starti],[nodes_N storelength],permute(store.ustar,[2 1]))
netcdf.putVar(ncid,nc_obj.var_waterSurf   ,[0 starti],[nodes_N storelength],permute(store.waterd,[2 1]))

% two dimensional var (nodesN, GSD)
netcdf.putVar(ncid,nc_obj.var_SurfFsi     ,[0 0 starti],[nodes_N gsd_MG+1 storelength],permute(store.SurfFsi,[2 3 1]))
netcdf.putVar(ncid,nc_obj.var_TranspPbi   ,[0 0 starti],[nodes_N gsd_MG storelength]  ,permute(store.TranspQbi,[2 3 1]))

if storei == 1
    % only save the first time
    netcdf.putVar(ncid,nc_obj.var_current_run       ,0,length(char(store.current_run)),store.current_run)
    netcdf.putVar(ncid,nc_obj.var_run_name          ,0,length(char(store.run_name)),store.run_name)
    % add HOST?
    % add EXECUTION PATH
    %netcdf.putVar(ncid,nc_obj.var_param_dt         ,0,1,store.param_dt)
    %netcdf.putVar(ncid,nc_obj.var_reach_length     ,0,1,store.param_reachlength)
    %netcdf.putVar(ncid,nc_obj.var_dx_array         ,0,nodes_N,store.dx_array)
    %netcdf.putVar(ncid,nc_obj.var_reach_width      ,0,nodes_N,store.reachwidth)
    %netcdf.putVar(ncid,nc_obj.var_gsd_initDsi      ,0,gsd_MG+1,store.gsd_initDsi)
end

if nc_params.saveNC_subs
    % subsurface vars:
    try
        % two dimensional var (nodesN, LStrMat)
        %netcdf.putVar(ncid,nc_obj.var_Dssg_subsurface     ,[0 0 starti],[LstrMat nodes_N storelength],permute(store.SubsDg,[2 3 1]))
        %netcdf.putVar(ncid,nc_obj.var_Dss90_subsurface    ,[0 0 starti],[LstrMat nodes_N storelength],permute(store.SubsD90,[2 3 1]))
        netcdf.putVar(ncid,nc_obj.var_eta_subs            ,[0 0 starti],[nodes_N LstrMat storelength],permute(store.SubsEta,[2 3 1]))
        
        % three dimensional var (nodesN, LstrMat, GSD)
        netcdf.putVar(ncid,nc_obj.var_pssi_subs    ,[0 0 0 starti],[nodes_N LstrMat gsd_MG storelength],permute(store.SubsPssi,[2 3 4 1]))
        
    catch err
        message = sprintf('ERROR while saving large NETCDF variables! Check cache size parameters, especially nc_params.premp');
        disp(message)
        error('Problem saving NetCDF')
    end
end

end

