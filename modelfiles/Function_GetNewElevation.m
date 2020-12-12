function [NodeGeom, NodeSedTrans] = Function_GetNewElevation( ...
    NodeGeom, NodeSedTrans, nodes_N, RunParam, timeOutput)
% Calculating new eta bed elevation, in m, function of the sediment

% Only add the feed of the first node into the upwind scheme. The rest of
% the feed should be added to the elevation directly.

% Reset arrays
sedTrans_array_dqbx = zeros(nodes_N,1);       %new Transport rate difference:
BedElev_change_Feed = zeros(nodes_N,1);
slope               = zeros(nodes_N,1);

slope_old           = NodeGeom.Slope;
sedTrans_array_qbx  = NodeSedTrans.qb;
feed_array_qin      = NodeSedTrans.qb_feed;
BedElev             = NodeGeom.Elev;
dx_array            = NodeGeom.dx;

param_au_upwind   = RunParam.Upwind_au;
param_au_downwind = 1 - param_au_upwind; %(1 - RunParam.Upwind_au);

% sedTrans_array_qbx(end) = 0;

switch RunParam.Upwind_order
    case 'first'
        %% FIRST ORDER UPWIND SCHEME https://en.wikipedia.org/wiki/Upwind_scheme
        debug = false;
        if debug
            format long g
            param_au_upwind   = 0.5;
            param_au_downwind = (1 - param_au_upwind); %(1 - RunParam.Upwind_au);
            sedTrans_array_qbx = [0;0;1;0;0];
            nodes_N             = length(sedTrans_array_qbx);
            sedTrans_array_dqbx = zeros(size(sedTrans_array_qbx));
            feed_array_qin      = zeros(size(sedTrans_array_qbx));
            dx_array            = ones(size(sedTrans_array_qbx));
        end
            
        % param_au should not affect feed from ghost node!
        sedTrans_array_dqbx(1) = ...
              ((param_au_upwind   * sedTrans_array_qbx(1)) - feed_array_qin(1))/dx_array(1)...
            + param_au_downwind * (sedTrans_array_qbx(2) - sedTrans_array_qbx(1))/dx_array(1);
        
        for n=2:nodes_N-1
            % upstream and downstream:
            U_xminus = param_au_upwind   * (sedTrans_array_qbx(n)   - sedTrans_array_qbx(n-1))/dx_array(n);
            U_xplus  = param_au_downwind * (sedTrans_array_qbx(n+1) - sedTrans_array_qbx(n))  /dx_array(n);
            sedTrans_array_dqbx(n) =  U_xminus + U_xplus;
            
        end
        
        % For this node first order:
        %sedTrans_array_dqbx(end - 1) = param_au_upwind * (sedTrans_array_qbx(end - 1) - sedTrans_array_qbx(end - 2))/dx_array(end - 1);
        %sedTrans_array_dqbx(end) = param_au_upwind * (sedTrans_array_qbx(end) - sedTrans_array_qbx(end - 1))/dx_array(end);
        sedTrans_array_dqbx(end) = (sedTrans_array_qbx(end) - sedTrans_array_qbx(end - 1))/dx_array(end);
        
        %sedTrans_array_dqbx(end) = 0;

    case 'second'
        % SECOND ORDER UPWIND SCHEME https://en.wikipedia.org/wiki/Upwind_scheme
        % first node (use first order, rest second exept last node):
        % param_au should not affect feed from ghost node!
        sedTrans_array_dqbx(1) = ...
              ((param_au_upwind   * sedTrans_array_qbx(1)) - feed_array_qin(1))/dx_array(1)...
            + param_au_downwind * (sedTrans_array_qbx(2) - sedTrans_array_qbx(1))/dx_array(1);
        
        sedTrans_array_dqbx(2) = ...
              param_au_upwind   * (  3 * sedTrans_array_qbx(2) - 4 * sedTrans_array_qbx(1) + 1 * feed_array_qin(1)    )/(2 * dx_array(2))...
            + param_au_downwind * (- 3 * sedTrans_array_qbx(2) + 4 * sedTrans_array_qbx(3) - 1 * sedTrans_array_qbx(4))/(2 * dx_array(2));
        
        for n=3:nodes_N-2
            % upstream and downstream:
            U_xminus = param_au_upwind   * (  3 * sedTrans_array_qbx(n) - 4 * sedTrans_array_qbx(n-1) + 1 * sedTrans_array_qbx(n-2))/(2 * dx_array(n));
            U_xplus  = param_au_downwind * (- 3 * sedTrans_array_qbx(n) + 4 * sedTrans_array_qbx(n+1) - 1 * sedTrans_array_qbx(n+2))/(2 * dx_array(n));
            sedTrans_array_dqbx(n) =  U_xminus + U_xplus;
        end
        
        % For this node first order:
        %if slope_old(end-1) < 0
        %    sedTrans_array_dqbx(end - 1) = (sedTrans_array_qbx(end - 1) - sedTrans_array_qbx(end - 2))/dx_array(end - 1);
        %else
%         sedTrans_array_dqbx(end - 1) = ...
%             param_au_upwind   * (sedTrans_array_qbx(end - 1) - sedTrans_array_qbx(end - 2))/dx_array(end - 1)...
%             + param_au_downwind * (sedTrans_array_qbx(end)     - sedTrans_array_qbx(end - 1))/dx_array(end - 1);
%         %end
%         % Fake ghost node: preserve differences for outflow conditions:
%         sedTrans_array_dqbx(end) = 0;
    otherwise
        error('unknown upwinding scheme specified!')
end

% check for mass conservation
%testsum = cumsum(sedTrans_array_dqbx);
%assert(testsum(end) < 1E-6,'assert diffusion working...')

% Elevation change due to Transport:
BedElev_change_Transport = (-sedTrans_array_dqbx./(1 - RunParam.Lambda)).*RunParam.dt_actual*RunParam.freqQ;

% Elevation change due to center feed, first node feed is already included in Transport:
BedElev_change_Feed(2:end) = (feed_array_qin(2:end)./dx_array(2:end)./(1 - RunParam.Lambda)).*RunParam.dt_actual*RunParam.freqQ;


% Calc new BedElev
BedElev_new = BedElev + BedElev_change_Transport + BedElev_change_Feed;
BedElev_new(end) = 0;
 
%% Calculation of bed slope, in m/m
switch RunParam.slopeOrder
    case 1
        % order one
        slope(1:end-1)= (BedElev_new(1:end-1) - BedElev_new(2:end))./(dx_array(1:end-1));
        slope(end)    = (BedElev_new(end-2)-BedElev_new(end))./(2 * dx_array(end));
    case 2
        % order two
        slope(1)      = (BedElev_new(1) - BedElev_new(2))./(dx_array(1));
        slope(2:end-1)= (BedElev_new(1:end-2) -BedElev_new(3:end))./(2*dx_array(2:end-1));
        slope(end)    = (BedElev_new(end-2)-BedElev_new(end))./(2 * dx_array(end));
    otherwise
        error('invalid SlopeOrder');
end
% 

% assign values to Output Structures:
% for n=1:nodes_N
NodeGeom.dElev_feed   = BedElev_change_Feed;

% and find first node impact:
% NodeGeom.dElev_feed(1) = ((feed_array_qin(1)/dx_array(1))./(1 - RunParam.Lambda)).*RunParam.dt_actual*RunParam.freqQ;

NodeGeom.Elev_new     = BedElev_new;
NodeGeom.Slope        = slope;
NodeSedTrans.dqb      = sedTrans_array_dqbx;
NodeSedTrans.dqb_feed = feed_array_qin./dx_array;
% end
%% CHECK MASS CONSERVATION
% feedIn          = feed_array_qin./dx_array./(1 - const_lambda).*param_dt_actual.*param_freqQ;
% transportChanges= sum(sedTrans_array_dqbx(1:end-1));
% %transportOut    = sedTrans_array_qbx(end)./dx_array(end)./(1 - const_lambda).*param_dt_actual.*param_freqQ;
% MassBefore      = sum(BedElev.*dx_array.*ChBankfullWidth_array) + sum(feedIn(2:end).*dx_array(2:end).*ChBankfullWidth_array(2:end)) + transportChanges.*dx_array(end)*ChBankfullWidth_array(end);
% MassAfter       = sum(BedElev_new.*dx_array.*ChBankfullWidth_array);
%
% assert(abs(1 - (MassBefore / MassAfter)) < 0.0001,'Mass conservation problem... in Function_GetNewElevationN ');


%% Error if problem
if or(isnan(sum(sedTrans_array_dqbx)),isnan(sum(BedElev_new)))
    errorposition   = 'Error:  t = %i, n = %i \n';
    errorcause      = 'isnan(sedTrans_array_dqbx) | isnan(BedElev_new)';
    %keyboard
    error(strcat(errorposition,errorcause), timeOutput,find(isnan(sedTrans_array_dqbx)))
end

end