function [ ncid, nc_obj,nc_params] = netcdf_create( current_run, NPrint,nodes_N,gsd_MG,LstrMat,nc_params )
%SAVE_ASNETC Summary of this function goes here
%   Detailed explanation goes here
nc_filename = strcat(current_run,'_store.nc');

%% SET NETCDF MODE
storenc_mode = netcdf.getConstant('NETCDF4');                           % 4096 - 1000000000000   % Create a NetCDF-4/HDF5 file
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('CLOBBER'));       %    0 - 0000000000000   % Overwrite any existing file with the same name.
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('NOCLOBBER'));     %    4 - 0000000000100   % Prevent overwriting of existing file with the same name.
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('CLASSIC_MODEL')); %  256 - 0000100000000   % Enforce the classic model; has no effect unless used in a bitwise-or with NETCDF4
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('64BIT_OFFSET'));  %  512 - 0001000000000   % Allow easier creation of files and variables which are larger than two gigabytes.
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('SHARE'));         % 2048 - 0100000000000   % Allow synchronous file updates.


%storenc_mode = netcdf.getConstant('64BIT_OFFSET');
%storenc_mode = bitor(storenc_mode,netcdf.getConstant('NETCDF4'));
% CREATE FILE
if exist(nc_filename,'file')
    try
        ncid = netcdf.open(nc_filename,'WRITE');
        netcdf.close(ncid);
        delete(nc_filename);
    catch err
        fprintf('Could not delete previous netcdf file....')
    end
end

try
    [ncid,nc_params.chunksize] = netcdf.create(nc_filename,storenc_mode);
catch err
    disp('NETCDF mode problem!')
    %keyboard
end
%% 
disp(strcat('Creating netcdf as: ',nc_filename));
nc_params.filename = nc_filename;

% testing....
% close netcdf
netcdf.close(ncid);
% open in write mode (so sync works??)
netcdf.open(nc_filename,'NC_WRITE')
% \testing!

if(nc_params.premp >= 1)
    message = sprintf('WARNING: nc_params.premp is >=1. I had problems with runs crashing while saving if nc_params.premp == 1. I set it to 1 for now... Expect crash...');
    fprintf(fdev,message);
    disp(message)
    warning(s'')
    nc_params.premp = min(1,nc_params.premp);
end
netcdf.setChunkCache(nc_params.csize, nc_params.nelems, nc_params.premp)
netcdf.setFill(ncid,nc_params.setFill);

%% SAVE DIMENSIONS
%ncid = netcdf.open('OAWC135000kg18000hGsigp50_7_store.nc','NC_WRITE');
netcdf.reDef(ncid)

nc_obj.dim_NPrint   = netcdf.defDim(ncid,'prints',NPrint+3);
nc_obj.dim_nodes    = netcdf.defDim(ncid,'nodes',nodes_N);
nc_obj.dim_gsd      = netcdf.defDim(ncid,'gsd',gsd_MG);
nc_obj.dim_gsdp1    = netcdf.defDim(ncid,'gsdp1',gsd_MG+1);
nc_obj.dim_const    = netcdf.defDim(ncid,'const',1);
nc_obj.dim_string   = netcdf.defDim(ncid,'string',200);
nc_obj.dim_LstrMat  = netcdf.defDim(ncid,'LstrMat',LstrMat);

% dimension vars (best practices)
nc_obj.var_NPrint   = netcdf.defVar(ncid,'prints','int',nc_obj.dim_NPrint);
nc_obj.var_nodes    = netcdf.defVar(ncid,'nodes','int',nc_obj.dim_nodes);
nc_obj.var_gsd      = netcdf.defVar(ncid,'gsd','int',nc_obj.dim_gsd);
nc_obj.var_gsdp1    = netcdf.defVar(ncid,'gsdp1','int',nc_obj.dim_gsdp1);
nc_obj.var_const    = netcdf.defVar(ncid,'const','int',nc_obj.dim_const);
nc_obj.var_string   = netcdf.defVar(ncid,'string','int',nc_obj.dim_string);
nc_obj.var_LstrMat  = netcdf.defVar(ncid,'LstrMat','int',nc_obj.dim_LstrMat);

%% CREATE VARIABLES
% zero dimensional var
nc_obj.var_t = netcdf.defVar(ncid,'store_t','int',nc_obj.dim_NPrint);
netcdf.putAtt(ncid,nc_obj.var_t,'unit','steps');
netcdf.putAtt(ncid,nc_obj.var_t,'description','model steps')
netcdf.defVarFill(ncid,nc_obj.var_t,false,NaN);
netcdf.defVarDeflate(ncid,nc_obj.var_t,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);

nc_obj.var_time = netcdf.defVar(ncid,'store_time','double',nc_obj.dim_NPrint);
netcdf.putAtt(ncid,nc_obj.var_time,'unit','hours');
netcdf.putAtt(ncid,nc_obj.var_time,'description','time in hours')
netcdf.defVarFill(ncid,nc_obj.var_time,false,NaN);
netcdf.defVarDeflate(ncid,nc_obj.var_time,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);

% one dimensional var (nodesN)
nc_obj.var_etab = netcdf.defVar(ncid,'store_etab','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_etab,'unit','m');
netcdf.defVarDeflate(ncid,nc_obj.var_etab,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_etab,false,NaN);

nc_obj.var_Dsg = netcdf.defVar(ncid,'store_Dsg','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_Dsg,'unit','mm');
netcdf.defVarDeflate(ncid,nc_obj.var_Dsg,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_Dsg,false,NaN);

nc_obj.var_Ds90 = netcdf.defVar(ncid,'store_Ds90','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_Ds90,'unit','mm');
netcdf.defVarDeflate(ncid,nc_obj.var_Ds90,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_Ds90,false,NaN);

nc_obj.var_SDg = netcdf.defVar(ncid,'store_SDg','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_SDg,'unit','stdev');
netcdf.defVarDeflate(ncid,nc_obj.var_SDg,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_SDg,false,NaN);

nc_obj.var_qbx = netcdf.defVar(ncid,'store_qbx','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_qbx,'unit','mps');
netcdf.defVarDeflate(ncid,nc_obj.var_qbx,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_qbx,false,NaN);

nc_obj.var_Msj = netcdf.defVar(ncid,'store_Msj','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_Msj,'unit','none');
netcdf.defVarDeflate(ncid,nc_obj.var_Msj,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_Msj,false,NaN);

nc_obj.var_slope = netcdf.defVar(ncid,'store_slope','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_slope,'unit','mpm');
netcdf.defVarDeflate(ncid,nc_obj.var_slope,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_slope,false,NaN);

nc_obj.var_ustar = netcdf.defVar(ncid,'store_ustar','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_ustar,'unit','none');
netcdf.defVarDeflate(ncid,nc_obj.var_ustar,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_ustar,false,NaN);

nc_obj.var_waterSurf = netcdf.defVar(ncid,'store_waterSurf','double',[nc_obj.dim_nodes nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_waterSurf,'unit','m');
netcdf.defVarDeflate(ncid,nc_obj.var_waterSurf,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_waterSurf,false,NaN);

% two dimensional var (nodesN, GSD)
nc_obj.var_TranspPbi = netcdf.defVar(ncid,'store_TranspPbi','double',[nc_obj.dim_nodes nc_obj.dim_gsd nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_TranspPbi,'unit','perc');
netcdf.defVarDeflate(ncid,nc_obj.var_TranspPbi,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_TranspPbi,false,NaN);

nc_obj.var_SurfFsi = netcdf.defVar(ncid,'store_SurfFsi','double',[nc_obj.dim_nodes nc_obj.dim_gsdp1 nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_SurfFsi,'unit','freq');
netcdf.defVarDeflate(ncid,nc_obj.var_SurfFsi,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_SurfFsi,false,NaN);

nc_obj.var_SurfPbi = netcdf.defVar(ncid,'store_SurfPbi','double',[nc_obj.dim_nodes nc_obj.dim_gsd nc_obj.dim_NPrint]);
netcdf.putAtt(ncid,nc_obj.var_SurfPbi,'unit','perc');
netcdf.defVarDeflate(ncid,nc_obj.var_SurfPbi,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
netcdf.defVarFill(ncid,nc_obj.var_SurfPbi,false,NaN);


% METADATA / DESCRIPTION OF NETCDF FILE
nc_obj.var_current_run = netcdf.defVar(ncid,'current_run','char',[nc_obj.dim_string]);
nc_obj.var_run_name = netcdf.defVar(ncid,'run_name','char',[nc_obj.dim_string]);
nc_obj.var_creation_date = netcdf.defVar(ncid,'creation_date','char',[nc_obj.dim_string]);
nc_obj.var_dx_array = netcdf.defVar(ncid,'dx_array','double',[nc_obj.dim_nodes]);
nc_obj.var_reach_width = netcdf.defVar(ncid,'reach_width','double',[nc_obj.dim_nodes]);

nc_obj.var_param_dt = netcdf.defVar(ncid,'param_dt','double',[nc_obj.dim_const]);
nc_obj.var_reach_length = netcdf.defVar(ncid,'reach_length','double',[nc_obj.dim_const]);
nc_obj.var_delta_s = netcdf.defVar(ncid,'param_delta_s','double',[nc_obj.dim_const]);
nc_obj.var_zbase = netcdf.defVar(ncid,'zbase','double',[nc_obj.dim_const]);

% inital values
nc_obj.var_gsd_initDsi = netcdf.defVar(ncid,'gsd_initDsi','double',[nc_obj.dim_gsdp1]);       
    
%% save subsurface if configured
if nc_params.saveNC_subs
    % subsurface (takes long time to store):
    % two dimensional var (nodesN, LStrMat)
    nc_obj.var_eta_subs = netcdf.defVar(ncid,'store_eta_subs','double',[nc_obj.dim_nodes nc_obj.dim_LstrMat nc_obj.dim_NPrint]);
    netcdf.putAtt(ncid,nc_obj.var_eta_subs,'unit','m');
    netcdf.defVarDeflate(ncid,nc_obj.var_eta_subs,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
    netcdf.defVarFill(ncid,nc_obj.var_eta_subs,false,NaN);
    
    % three dimensional var (nodesN, LstrMat, GSD)
    %     nc_obj.var_Dssg_subsurface = netcdf.defVar(ncid,'store_Dssg_subsurface','double',[nc_obj.dim_LstrMat nc_obj.dim_nodes nc_obj.dim_NPrint]);
    %     netcdf.putAtt(ncid,nc_obj.var_Dssg_subsurface,'unit','mm');
    %     netcdf.defVarDeflate(ncid,nc_obj.var_Dssg_subsurface,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
    %     netcdf.defVarFill(ncid,nc_obj.var_Dssg_subsurface,false,NaN);
    %
    %     nc_obj.var_Dss90_subsurface = netcdf.defVar(ncid,'store_Dss90_subsurface','double',[nc_obj.dim_LstrMat nc_obj.dim_nodes nc_obj.dim_NPrint]);
    %     netcdf.putAtt(ncid,nc_obj.var_Dss90_subsurface,'unit','mm');
    %     netcdf.defVarDeflate(ncid,nc_obj.var_Dss90_subsurface,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
    %     netcdf.defVarFill(ncid,nc_obj.var_Dss90_subsurface,false,NaN);
    
    nc_obj.var_pssi_subs = netcdf.defVar(ncid,'store_pssi_subs','double',[nc_obj.dim_nodes nc_obj.dim_LstrMat nc_obj.dim_gsd nc_obj.dim_NPrint]);
    netcdf.putAtt(ncid,nc_obj.var_pssi_subs,'unit','perc');
    netcdf.defVarDeflate(ncid,nc_obj.var_pssi_subs,nc_params.shufflefilter,nc_params.deflate,nc_params.deflate_level);
    netcdf.defVarFill(ncid,nc_obj.var_pssi_subs,false,NaN);
        
end

% DISPLAY CONTENT OF NETCDF FILE
%ncdisp(storenc_filename)

% Take file out of define mode.
netcdf.endDef(ncid);

%% SAVE DIMENSIONS (as per best practices)
netcdf.putVar(ncid,nc_obj.var_NPrint   ,0, NPrint+3 ,[1:NPrint+3])
netcdf.putVar(ncid,nc_obj.var_nodes    ,0, nodes_N  ,[1:nodes_N])
netcdf.putVar(ncid,nc_obj.var_gsd      ,0, gsd_MG   ,[1:gsd_MG])
netcdf.putVar(ncid,nc_obj.var_gsdp1    ,0, gsd_MG+1 ,[1:gsd_MG+1])
netcdf.putVar(ncid,nc_obj.var_const    ,0, 1        ,1)
netcdf.putVar(ncid,nc_obj.var_string   ,0, 200      ,[1:200])
netcdf.putVar(ncid,nc_obj.var_LstrMat  ,0, LstrMat  ,[1:LstrMat])
    
end

