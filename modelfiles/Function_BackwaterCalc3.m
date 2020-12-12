function [ustar_array,waterDepth,froudeArray, increaseSubtimestep] = Function_BackwaterCalc3(qw,r,lowerReachFlowDepth,...
    NodeGeom, NodeGSD, RunParam, nodes_id, increaseSubtimestep)
format long g


%% fake input to understand the script without calling it from the model:
% debug = false;
% if debug
%     
%     qw                      = 0.0650;
%     % flat downslope:
%     slope_array             = ones(41,1) .* 0.05;
%     slope_array(10:18)      = -0.01;
%     
%     % bumpy:
%     slope_array             = transpose([0.01 0.01 0.01 0.01 0.01 0.05 0.05 0.05 ...
%         0.05 0.05 0.07 0.07 0.05 0.03 0.02 0.02 ...
%         0.01 -0.01 -0.02 -0.01 0.01 0.01 0.03 0.03 ...
%         0.03 0.03 0.03 0.03 0.03 0.03 0.07 0.07 ...
%         0.01 0.01 0.01 0.01 0.001 0.001 0.001 0.001 ...
%         0.001]);
%     
%     gsd_arraySurfD90_new    = ones(41,1) .* 23.425;
%     dx_array                = ones(41,1) .* 0.5;
%     nodes_loc               = zeros(length(dx_array),1);
%     nodes_loc(2:end)        = cumsum(dx_array(1:end-1));
%     param_nk                = 2;
%     param_alfa_r            = 8.1;
%     const_rhoW              = 1000;
%     const_g                 = 9.81;
%     
%     NodeGeom.Elev                 = NaN(length(dx_array),1);
%     NodeGeom.Elev(end) = 0; % outflow
%     for n = length(NodeGeom.Elev)-1:-1:1
%         NodeGeom.Elev(n) = NodeGeom.Elev(n+1) + tan(slope_array(n))*dx_array(n);
%     end
%     
%     backwaterparam.iter_max         = 1000;     % maximum number of iterations to reach iter_diffgood:
%     backwaterparam.iter_diffgood    = 0.00001;   % stop iterating waterdepth approximation if difference between calculations is smaller than this value
%     backwaterparam.waterdepth_min   = 0.00001;   % calc finer spatial grid if waterdepth is calculated to be smaller than value. maybe calculate this threshold realistically?
%     backwaterparam.froudefact_min   = 0.8;      % calc finer spatial grid if upstream decrease of froudenumber in subcritical condition is calculated to be smaller than this value multiplied by downstream Froude number.
%     backwaterparam.subGridScalemax  = 65536;    % maximum interpolation factor for calculating subgrids. Error if tripped
%     backwaterparam.Fr_threshold     = 0.9;      % Froude Number Threshold for using supercritical:
%     backwaterparam.addGhostNode     = false;
%     nodes_id = [1:length(slope_array)];
%     
%     NodeGeom.Slope = slope_array;
%     NodeGSD.SurfD90_new= gsd_arraySurfD90_new;
%     NodeGeom.dx = dx_array;
%     NodeGeom.WidthFixed = ones(size(slope_array));
% end

debug = false;
if debug
    load('matlab.mat');
    debug = true;
    
    
    %     NodeGeom.Elev  = rundata.NodeGeom.Elev(204:233);
    %     qw             = 0.3932;
    %     NodeGeom.Slope = rundata.NodeGeom.Slope(204:233);
    %     NodeGeom.dx    = rundata.NodeGeom.dx(204:233);
    %     NodeGSD.SurfD90_new  = rundata.NodeGSD.SurfD90_new(204:233);
    %     nodes_id = [1:length(NodeGeom.Slope)];
    %     NodeGeom.WidthFixed = rundata.NodeGeom.WidthFixed(204:233);
end
%%

param_nk        = RunParam.nk;
param_alfa_r    = RunParam.alfar;
const_rhoW      = RunParam.RhoW;
const_g     	= RunParam.Gravity;
backwaterparam  = RunParam.backwaterparam;

if ~isfield(backwaterparam,'addGhostNode')
    % default to false if unspecified:
    backwaterparam.addGhostNode = false;
end
if ~isfield(backwaterparam,'bisect_max')
    backwaterparam.bisect_max = 10E6;
end

lowerReachExists = (nodes_id(end) < size(NodeGeom.Slope,1));

addGhostNode = backwaterparam.addGhostNode && lowerReachExists;

if addGhostNode
    % add a ghost node at the end:
    slope_array = [NodeGeom.Slope(nodes_id);  NodeGeom.Slope(nodes_id(end) + 1)];
    SurfD90     = [NodeGSD.SurfD90_new(nodes_id); NodeGSD.SurfD90_new(nodes_id(end) + 1)];
    dx_array 	= [NodeGeom.dx(nodes_id);     NodeGeom.dx(nodes_id(end) + 1)];
    if slope_array(end) < 10E-4
        slope_array(end) = 10E-4;
    end
else
    slope_array = [NodeGeom.Slope(nodes_id)];
    SurfD90     = [NodeGSD.SurfD90_new(nodes_id)];
    dx_array 	= [NodeGeom.dx(nodes_id)];
end

% Backwater curve calcutation
% Nikuradse coefficient, in m
ks = param_nk.*SurfD90(:)./1000;

% Set up arrays
nodes_N             = length(slope_array);
%taub                = NaN(size(slope_array));
frictionSlopeHere   = NaN(size(slope_array));
waterDepth          = NaN(size(slope_array));
froudeArray         = NaN(size(slope_array));

% Critical water depth, in m
waterDepth_critical = (qw^2/const_g)^(1/3);
waterDepth_critical2  = (qw/const_g^(1/2) ./ backwaterparam.Fr_threshold).^ (2/3);

% Critical slope:
slope_critical = (ks(end)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_critical2^(10/3));

% Assert that the last node does not have very high water depths
if slope_array(end) < 10E-4
    slope_array(end) = 10E-4;
end

% Normal Water depth Calculation (Manning-strickler), in m
waterDepth_normal=((ks(end)^(1/3)*qw^2)/(param_alfa_r^2*const_g*slope_array(end)))^(3/10);

% We start calculation at the outlet (downstream boundary condition)
if RunParam.LagoonStage_fixed && any(RunParam.lagoonNode == nodes_id)
    waterDepth_downstreamBoundary = RunParam.LagoonStage_ft * RunParam.CFStoCMS^(1/3);
    %waterDepth_downstreamBoundary = max(waterDepth_normal,waterDepth_downstreamBoundary);
else
    waterDepth_downstreamBoundary = waterDepth_normal;
end

if lowerReachExists
    % Get waterdepth from downstream:
    waterDepth(end) = lowerReachFlowDepth;
else
    % DOWNSTREAM END if simulation
    if waterDepth_downstreamBoundary >= (waterDepth_normal + 0.0001)
        waterDepth(end) = waterDepth_downstreamBoundary;
    else
        %waterDepth(end) = waterDepth_critical2 + 0.001; % carles has it a bit higher +0.1;
        waterDepth(end) = waterDepth_normal + 0.0001;
    end
end
frictionSlopeHere(end) = (ks(end)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth(end)^(10/3));
froudeArray(end) = qw/const_g^(1/2)/waterDepth(end)^(3/2);

fillPoolShawn = true;

slope_array_old = slope_array;

if fillPoolShawn
    fillPoolCount = 0;
    while any(slope_array<=0)
        %elev_test  = NaN(length(dx_array),1);
        %elev_test(end) = 0; % outflow
        %for n = length(elev_test)-1:-1:1
         %   elev_test(n) = elev_test(n+1) + tan(slope_array(n))*dx_array(n);
        %end
        %elev_test = 0;
        
        % Do Shawns Pool-preserving slope adjustment
        % identify pool and find new slope for water surface slope calculation
        poolNode = find(slope_array<=0,1,'last');
        
        % find downstream end of pool node:
        node_poolend_reachid    = poolNode-1 + find(slope_array(poolNode:end) > 0,1,'first') ;
        node_poolend_allid = node_poolend_reachid + nodes_id(1) - 1;
        
        % find upstream beginning of pool node:
        node_poolstart_allid  = find(NodeGeom.Elev(1:node_poolend_allid-1) > NodeGeom.Elev(node_poolend_allid),1,'last');
        if ~isempty(node_poolstart_allid)
            if(node_poolstart_allid < nodes_id(1))
                node_poolstart_reachid = 1;
            else
                node_poolstart_reachid = find(nodes_id == node_poolstart_allid,1,'first');
            end
            slopeOverPool = (NodeGeom.Elev(node_poolstart_allid) - NodeGeom.Elev(node_poolend_allid))/sum(NodeGeom.dx(node_poolstart_allid:node_poolend_allid-1));
        elseif slope_array(1)<0
            node_poolstart_reachid = 1;
            slopeOverPool = 10E-4;
        end
        
        
        assert(slopeOverPool > 0,'Pool water surface slope not greater 0!');
        
        %slope_array(firstNodeHigher + 1:poolNode) = slopeToHigherNode;
        slope_array(node_poolstart_reachid:poolNode) = slopeOverPool;
        
        %warning(['Found 0 or smaller slope at: ',num2str(poolNode(p))])
        
        % error condition:
        fillPoolCount = fillPoolCount + 1;
        assert(fillPoolCount < 100,'Filling pools does not work...');
    end
    assert(all(slope_array>0),'Not all slope values are gt 0!')
end

%
redoWithSubStep = false;

% Choose a stepping factor of 1 channel width minimum!
%bisect_fact_start = ceil(max(dx_array(nodes_N))/max(NodeGeom.WidthFixed(nodes_id)));
%bisect_fact_start = max(bisect_fact_start,10);

bisect_fact_start = ones(length(dx_array),1) .* 10;
%bisect_fact_start(end-1:end) = deal(2); % boundary problem...
bisect_fact_start(end) = deal(2); % boundary problem...

bisect_fact = bisect_fact_start;

dx_array_bisect = dx_array(nodes_N);

waterDepth_subhere      = NaN;
frictionSlope_subhere   = NaN;
froude_subhere          = NaN;
slope_array_subgrid    = [slope_array(nodes_N), slope_array(nodes_N - 1)];

bisect_max = backwaterparam.bisect_max;
n = nodes_N - 1;
%bisect_max = 200;

%backwaterparam.iter_max = 100;
%backwaterparam.bisect_max = 10000;

% initial values
% Calculate initial approximation of Froude number by using
% downstream waterdepth / boundary waterdepth for downstream end:
% substep parameters
froude_subbefore        = qw/const_g^(1/2)/waterDepth(n+1)^(3/2);
waterDepth_subbefore    = waterDepth(n+1);
frictionSlope_subbefore = frictionSlopeHere(n+1);
stopReason = 'none';

stopReasonSave = [];
stopReasonSaveNo = 1;
%%
while n > 0
    % iterate through nodes backwards, from downstream to upstream
    if debug
        disp([num2str(n),' Calculating Backw, FroudeBefore: ',num2str(froude_subbefore)]);
        if redoWithSubStep
            disp([num2str(n),' Reason for redo: ',stopReason]);
        end
    end
    iter = 1;
    stopLoop2 = false;
    
    mean_froude = 0;
    %     mean_water  = 0;
    %     mean_fslope = 0;
    
    % temporarily save parameters in case we have to redo step:
    if ~redoWithSubStep
        store_froude = froude_subbefore;
        store_water  = waterDepth_subbefore;
        store_fslope = frictionSlope_subbefore;
        store_froude_subhere = froude_subhere;
    else
        % Bisection: generate more nodes between the downstream and this
        % node...
        bisect_fact(n) = bisect_fact(n) * 2;
        
        stopReasonSave{stopReasonSaveNo} = stopReason;
        stopReasonSaveNo = stopReasonSaveNo + 1;
        if bisect_fact(n) > bisect_max
            % remove big variables from workspace, so we can take a
            % snapshot
            
            clear slope_array_subgrid slopeinterp distanceDx
            save('backw_error.mat')
            
            %% original-node profile:
            figure;plot(NodeGeom.Elev(nodes_id) - NodeGeom.Elev(nodes_id(end)),'-*k');
            hold on;
            plot(NodeGeom.Elev(nodes_id)+waterDepth(1:length(nodes_id)) - NodeGeom.Elev(nodes_id(end)),'-*b');
            hold on;
            froudeSuperCrit = froudeArray .* (froudeArray > backwaterparam.Fr_threshold);
            froudeSuperCrit(froudeSuperCrit == 0) = deal(NaN);
            froudeSubCrit   = froudeArray .* (froudeArray < backwaterparam.Fr_threshold);
            froudeSubCrit(froudeSubCrit == 0) = deal(NaN);
            yyaxis right
            plot(froudeSuperCrit,'-*r')
            hold on;
            plot(froudeSubCrit,'-or')
            hold on;
            plot(slope_array./max(slope_array),'-m')
%             hold on;
%             hplot(ustar_array./max(ustar_array),'-y')
            
            %% throw error
            error(['bisect_fact is running away', stopReason])
        end
        
        % load last valid values
        froude_subbefore = store_froude;
        waterDepth_subbefore = store_water;
        frictionSlope_subbefore = store_fslope;
        redoWithSubStep = false;
        iter = 1;
    end
    
    if (slope_array(n) <= 0 )
        waterDepth_subhere = waterDepth(n+1) - slope_array(n) * dx_array(n);
        %waterDepth_subhere = waterDepth_critical2;
        frictionSlope_subhere  = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subhere^(10/3));
        froude_subhere         = qw/const_g^(1/2)./waterDepth_subhere.^(3/2);
        
        mean_froude = froude_subhere;
        bisect_fact(n) = 2;
        
        waterDepth_subbefore     = waterDepth_subhere;
        frictionSlope_subbefore  = frictionSlope_subhere;
        froude_subbefore         = froude_subhere;
        
        redoWithSubStep = false;
    elseif (slope_array(n+1) <= 0 ) && (slope_array(n) >= slope_critical)
        waterDepth_subhere     = ((ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*slope_array(n)))^(3/10);
        frictionSlope_subhere  = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subhere^(10/3));
        froude_subhere         = qw/const_g^(1/2)./waterDepth_subhere.^(3/2);
        
    else
        % get distance for this calculation
        dx_array_bisect       = dx_array(n) / bisect_fact(n);
        
        
        slopeinterp = ones(1,bisect_fact(n)) .*slope_array(n);
        if bisect_fact(n) > 1
            % interpolate
            slopeOfslope = slope_array(n+1) - slope_array(n);
            distanceDx   = flip([dx_array_bisect:dx_array_bisect:dx_array(n)-dx_array_bisect]./dx_array(n));
            slopeinterp(1:end-1) = slopeOfslope .* distanceDx + slope_array(n);
        end
        slope_array_subgrid  = horzcat(slope_array(n+1), slopeinterp);
        
        % old build in interpolation is hundreds times slower
        %slope_array_subgrid2   = interp1([0,dx_array(n)], [slope_array(n+1), slope_array(n)],[0:dx_array_bisect:dx_array(n)]);
        %assert(any(abs(slope_array_subgrid-slope_array_subgrid2)< 10E-6),'different!')
        
        substeps = 0;
        nn=2; % start from second point, first one is the previous n
        while nn <= length(slope_array_subgrid)
            if ( froude_subbefore >= backwaterparam.Fr_threshold) && (slope_array_subgrid(nn) > 0)
                % SUPERCRITICAL flow: normal flow approximation
                % approximation for negative slopes:
                
                
                waterDepth_subhere     = ((ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*slope_array_subgrid(nn)))^(3/10);
                frictionSlope_subhere  = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subhere^(10/3));
                froude_subhere         = qw/const_g^(1/2)./waterDepth_subhere.^(3/2);
                
                if ~isreal(froude_subhere)
                    redoWithSubStep = true;
                    stopReason = '(1) froude nonreal in supercritical flow';
                    break
                elseif froude_subhere < froude_subbefore * backwaterparam.froudefact_min
                    % FROUDE CONDITION:
                    % Allowed: smooth upstream decrease of Froude: supercritical --> critical --> subcritical
                    % Allowed: abrupt upstream increase of Froude: subcritical --> hydr.jump --> supercritical
                    % Problem D) abrupt upstream decrease of Froude: supercritical --> subcritical
                    if n ~= nodes_N && (abs(froude_subhere - store_froude_subhere) > 0.1)
                        % redo with higher resolution
                        redoWithSubStep = true;
                        stopReason = ['(5a) too steep decrease in Fr between substep and before! from ',num2str(froude_subbefore),' to ',num2str(froude_subhere)];
                        store_froude_subhere = froude_subhere;
                        break
                    else
                        % iteration did not improve the value!
                        stopLoop2  = true;
                        redoWithSubStep = false;
                        froude_subbefore         = froude_subhere;
                        break
                    end
                else
                    redoWithSubStep     = false;
                    waterDepth_subbefore     = waterDepth_subhere;
                    frictionSlope_subbefore  = frictionSlope_subhere;
                    froude_subbefore         = froude_subhere;
                    nn = nn + 1;
                end
                
                if debug
                    disp([num2str(n),'++supercrit flow nn: ',num2str(nn),'  Fr=',num2str(froude_subhere)])
                end
            else
                if ( froude_subbefore > backwaterparam.Fr_threshold) && (slope_array_subgrid(nn) <= 0)
                    %froude_subbefore = backwaterparam.Fr_threshold - 0.1;
                    waterDepth_subbefore = waterDepth_critical2;
                    %frictionSlope_subbefore  = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subbefore^(10/3));
                    frictionSlope_subbefore  = slope_critical;
                    froude_subbefore         = backwaterparam.Fr_threshold;
                end
                % SUBCRITICAL flow
                froude_subhere = froude_subbefore;
                frictionSlope_subhere = frictionSlope_subbefore;
                
                % Calculate the water hight at this point nn:
                while ~stopLoop2
                    dHdX_i = (slope_array_subgrid(nn-1) - frictionSlope_subhere)/(1 - froude_subhere^2);
                    waterDepth_subhere = waterDepth_subbefore - dHdX_i*dx_array_bisect;
                    frictionSlope_subhere = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subhere^(10/3));
                    froude_subhere = qw/const_g^(1/2)/waterDepth_subhere^(3/2);
                    
                    if ~isreal(waterDepth_subhere * dHdX_i * frictionSlope_subhere * froude_subhere)
                        % Problem B) any of waterDepth_subhere, dHdX_i, frictionSlope_subhere, froude_subhere is nonreal/imaginary number (happens if negative number^(10/3))
                        % redo with higher resolution
                        stopLoop2  = true;
                        redoWithSubStep = true;
                        stopReason = '(2) variable nonreal in subcritical flow solver';
                        froude_subbefore         = froude_subhere;
                        break
                    end
                    
                    % check if something went out of bound:
                    if iter > backwaterparam.iter_max
                        % Problem A) too many iterations to find a value within iter_diffgood between two successive iterations
                        % redo with higher resolution
                        stopLoop2  = true;
                        redoWithSubStep = true;
                        stopReason = ' (3) too many iterations in subcritical flow solver';
                        break
                    end
                    
                    % check how much the additional iteration changed the
                    % value, if below threshold, then keep:
                    if iter>2 && abs(dHdX_i - dHdX_old) < backwaterparam.iter_diffgood
                        % found value, go to next node
                        stopLoop2  = true;
                        redoWithSubStep = false;
                        break
                    end
                    % iterate
                    iter = iter + 1;
                    dHdX_old = dHdX_i;
                end
                %             end
                
                % Find maximum water depth in case of hydraulic jump
                % location
                if (slope_array(n+1) < 0) && slope_array(n) > slope_critical
                    maxWaterdepth_super = ((ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g* slope_array(n)))^(3/10);
                    waterDepth_subhere = min(waterDepth_subhere,maxWaterdepth_super);
                    frictionSlope_subhere = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_subhere^(10/3));
                    froude_subhere = qw/const_g^(1/2)/waterDepth_subhere^(3/2);
                    
                    mean_froude = mean_froude + froude_subhere;
                    %bisect_fact(n) = 1;
                    break
                end
                
                if debug
                    disp([num2str(n),'----subcrit flow nn: ',num2str(nn),'  Fr=',num2str(froude_subhere)])
                end
                
                % reset after loop
                stopLoop2  = false;
                
                if ~redoWithSubStep
                    % check if this data is valid
                    %             if waterDepth_subhere < backwaterparam.waterdepth_min
                    %                 % Problem C) waterDepth_subhere is very small / negative
                    %                 redoWithSubStep = true;
                    %                 stopLoop1 = false;
                    %             elseif froude_subhere < Froude_DS_Thresh
                    
                    %                     if froude_subhere > 1
                    %                         % Problem: Froude number should be subcritical, not
                    %                         % supercritical! Redo with finer resolution
                    %                         redoWithSubStep = true;
                    %                         stopReason = ' (4) switching from subcritical to supercritical in subnode';
                    %                         break
                    %                     else
                    if froude_subhere < froude_subbefore * backwaterparam.froudefact_min
                        % FROUDE CONDITION:
                        % Allowed: smooth upstream decrease of Froude: supercritical --> critical --> subcritical
                        % Allowed: abrupt upstream increase of Froude: subcritical --> hydr.jump --> supercritical
                        % Problem D) abrupt upstream decrease of Froude: supercritical --> subcritical
                        if isnan(store_froude_subhere)
                            store_froude_subhere = froude_subbefore;
                        end
                        if (n ~= nodes_N) && (abs(froude_subhere - store_froude_subhere) > 0.01)
                            % redo with higher resolution
                            redoWithSubStep = true;
                            stopReason = ' (5b) too steep decrease in Fr between substep and before!';
                            store_froude_subhere = froude_subhere;
                            break
                        else
                            store_froude_subhere = froude_subhere;
                        end
                    elseif froude_subhere > froude_subbefore * 2
                        redoWithSubStep = true;
                        stopReason = ' (5c) too steep increase in Fr between substep and before!';
                        if froude_subhere > backwaterparam.Fr_threshold
                           % switch to supercritical 
                           stopReason = ' (5d) switching to supercritical!';
                           store_froude = backwaterparam.Fr_threshold;
                        end
                        store_froude_subhere = froude_subhere;
                        break
                    else
                        % store this value for the subgrid
                        waterDepth_subbefore     = waterDepth_subhere;
                        frictionSlope_subbefore  = frictionSlope_subhere;
                        froude_subbefore         = froude_subhere;
                        nn = nn + 1;
                    end
                else
                    
                    if debug
                        disp([num2str(n),' BREAK subcrit flow nn: ',num2str(nn),'  Fr=',num2str(froude_subhere)])
                    end
                    break
                end
            end
            mean_froude = mean_froude + froude_subhere;
            substeps = substeps + 1;
        end
    end
    if ~redoWithSubStep
        % Subspatial grid done
        %if slope_array(n) < 0
        %    waterDepth(n) = waterDepth_subhere - slope_array(n)*dx_array(n);
        %else
        waterDepth(n) = waterDepth_subhere;
        
        %end
        frictionSlopeHere(n) = frictionSlope_subhere;
        
        froudeArray(n)       = mean_froude / bisect_fact(n);
        froudeArray(n)       = froude_subhere;
        
        %froude_subbefore         = mean_froude;
        %         waterDepth_subbefore     = waterDepth_subhere;
        %         frictionSlope_subbefore  = mean_fslope;
        
        % prepare next time step
        n = n - 1;
        if n == 0
            break
        end
        redoWithSubStep = false;
        bisect_fact     = bisect_fact_start;
        
        %dx_array_bisect = dx_array(n) / bisect_fact;
        %slope_array_subgrid    = [slope_array(n+1), slope_array(n)];
    end
    
    %assert(isreal(froudeArray(n)) && froudeArray(n) > 0 ,'froudenumber too small / zero')
end

%
if ~isreal(waterDepth)
    error('waterDepth gave non-real number!')
end

% incrWD = waterDepth(1:end-1)./waterDepth(2:end);
% incrWD2= incrWD .* (slope_array(2:end) > 0);
% incrWD3= incrWD .* ((slope_array(1:end-1) > 0) .* (froudeArray(1:end-1) > 0));

%assert(all(incrWD < 3),'Issue: Water depth change too large!')

% Calculate resulting shear stress
taub = const_rhoW .* const_g .* waterDepth .* frictionSlopeHere; % Shear stress

% set to zero in nodes with negative bedslope:
ustar_array = (slope_array_old > 0) .* sqrt(taub./const_rhoW); % Shear velocity

%ustar_array = sqrt(taub./const_rhoW); % Shear velocity


% if any(ustar_array>1)
%     warning('Extreme shear velocity!')
%     %keyboard
% end

% don't set to zero in nodes with negative bedslope:
%ustar_array = sqrt(taub./const_rhoW); % Shear velocity

%if(~all(incrWD2 < 10))
if debug
    %% original-node profile:
    figure;plot(NodeGeom.Elev(nodes_id) - NodeGeom.Elev(nodes_id(end)),'-*k');
    hold on;
    plot(NodeGeom.Elev(nodes_id)+waterDepth(1:length(nodes_id)) - NodeGeom.Elev(nodes_id(end)),'-*b');
    hold on;
    froudeSuperCrit = froudeArray .* (froudeArray > backwaterparam.Fr_threshold);
    froudeSuperCrit(froudeSuperCrit == 0) = deal(NaN);
    froudeSubCrit   = froudeArray .* (froudeArray < backwaterparam.Fr_threshold);
    froudeSubCrit(froudeSubCrit == 0) = deal(NaN);
    yyaxis right
    plot(froudeSuperCrit,'-*r')
    hold on;
    plot(froudeSubCrit,'-or')
    hold on;
    plot(slope_array./max(slope_array),'-m')
    hold on;
    plot(ustar_array./max(ustar_array),'-y')
    %%
    %keyboard
end


%% debug
%keyboard
% 
% if(~all(incrWD2 < 2))
%     warning('Issue: Water depth change too large!')
%     keyboard
% end

if addGhostNode && ~debug
    % remove ghost node for reporting of values:
    ustar_array   = ustar_array(1:end-1);
    waterDepth    = waterDepth(1:end-1);
    froudeArray  = froudeArray(1:end-1);
end


if any(isnan(ustar_array))
    error('shear velocity error')
end

if any(ustar_array<0)
    error('ustar is negative!')
end

if any(ustar_array > 2)
    warningNode = find(max(ustar_array)==ustar_array,1);
    if nodes_id(warningNode) > 2
        %debug = true;
        warning(['Node ',num2str(nodes_id(warningNode)),': extreme ustar in backw ',...
            'Froude: ',num2str(froudeArray(warningNode))])
    end
    increaseSubtimestep = true;
    %keyboard
end

if debug
    %% original-node profile:
    figure;plot(nodes_id,NodeGeom.Elev(nodes_id) - NodeGeom.Elev(nodes_id(end)),'-*k');
    hold on;
    plot(nodes_id,NodeGeom.Elev(nodes_id)+waterDepth(1:length(nodes_id)) - NodeGeom.Elev(nodes_id(end)),'-*b');
    hold on;
    froudeSuperCrit = froudeArray .* (froudeArray > backwaterparam.Fr_threshold);
    froudeSuperCrit(froudeSuperCrit == 0) = deal(NaN);
    froudeSubCrit   = froudeArray .* (froudeArray < backwaterparam.Fr_threshold);
    froudeSubCrit(froudeSubCrit == 0) = deal(NaN);
    yyaxis right
    if any(~isnan(froudeSuperCrit))
        plot(nodes_id,froudeSuperCrit,'-*r')
        hold on;
    end
    if any(~isnan(froudeSubCrit))
        plot(nodes_id,froudeSubCrit,'-or')
        hold on;
    end
    plot(nodes_id,NodeGeom.Slope(nodes_id)./max(NodeGeom.Slope(nodes_id)),'-m')
    hold on;
    plot(nodes_id,ustar_array./max(ustar_array),'-y')
    %%
    keyboard
end
end