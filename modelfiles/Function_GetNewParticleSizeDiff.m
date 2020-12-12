function [ NodeGSD,NodeSedTrans ] = Function_GetNewParticleSizeDiff( NodeGeom, NodeGSD, NodeSedTrans, NodeSubs,RunParam,t_step )
%FUNCTION_GETNEWPARTICLESIZEDIFF Summary of this function goes here
%   Detailed explanation goes here

% update frequencies of first subsurface layer:
pssi     = squeeze(NodeSubs.Pfi(:,2,:));
nodes_N  = length(RunParam.nodes_loc);
gsd_MG   = size(NodeGSD.SurfFsi,2)-1;
dx_array = RunParam.dx_array;

% initialize arrays:
fIi             = zeros(size(pssi));
sedTrans_qbxi   = zeros(size(NodeSedTrans.qbi));
feed_qbxi       = zeros(size(NodeSedTrans.qbi));

%% Calculating the exchange frequencies between the active layer and the substrate
LogicArr_Degr = [NodeGeom.Elev_new] < [NodeGeom.Elev];
LogicArr_Aggr = ~LogicArr_Degr;
LogicArr_Aggr_Feed = logical(repmat(LogicArr_Aggr .* ([NodeGeom.dElev_feed] > 0),1,gsd_MG));
LogicArr_Degr_Feed = logical(repmat(LogicArr_Degr .* ([NodeGeom.dElev_feed] > 0),1,gsd_MG));
LogicArr_Aggr_NoFeed = logical(repmat(LogicArr_Aggr .* ([NodeGeom.dElev_feed] == 0),1,gsd_MG));
LogicArr_Degr_NoFeed = logical(repmat(LogicArr_Degr .* ([NodeGeom.dElev_feed] == 0),1,gsd_MG));
assert(sum(LogicArr_Aggr_Feed(:))+sum(LogicArr_Degr_Feed(:))+sum(LogicArr_Aggr_NoFeed(:))+sum(LogicArr_Degr_NoFeed(:))==nodes_N*gsd_MG,'Does not sum up with node number!')
% Degrading conditions
fIi(LogicArr_Degr_NoFeed) = pssi(LogicArr_Degr_NoFeed); % take frequencies from first subsurface layer

% add a little bit randomness to prevent crashes from weird pfi:
alphaF_plusRandom = RunParam.alfaF + (rand(1) / 100);

% Aggrading or stagnant conditions:
fIi(LogicArr_Aggr_NoFeed) = alphaF_plusRandom.*NodeGSD.SurfPfi(LogicArr_Aggr_NoFeed) ...
    + (1 - alphaF_plusRandom)*NodeSedTrans.qbi(LogicArr_Aggr_NoFeed);

% # Added Pulse feed impact here!!!! #
% Both Transport and PulseFeed aggrade:
PulseImpactRatio_Aggr = repmat((([NodeGeom.Elev_new]-[NodeGeom.Elev])./[NodeGeom.dElev_feed]),1,gsd_MG);
PulseImpactRatio_Aggr(isnan(PulseImpactRatio_Aggr))   = 0; % no aggradation comes from feed
PulseImpactRatio_Aggr(isinf(PulseImpactRatio_Aggr))   = 0;
PulseImpactRatio_Aggr(PulseImpactRatio_Aggr>1)        = 1; % all aggradation comes from feed
PulseImpactRatio_Aggr(PulseImpactRatio_Aggr<0)        = 0; % all aggradation comes from feed

%         PulseImpactRatio_Degr = repmat((([NodeGeom.Elev]-[NodeGeom.Elev_new])./[NodeGeom.dElev_feed]),1,gsd_MG);
%         PulseImpactRatio_Degr(isnan(PulseImpactRatio_Degr))   = 0; % no aggradation comes from feed
%         PulseImpactRatio_Degr(isinf(PulseImpactRatio_Degr))   = 0;
%         PulseImpactRatio_Degr(PulseImpactRatio_Degr>1)        = 1; % all aggradation comes from feed
%         PulseImpactRatio_Degr(PulseImpactRatio_Degr<0)        = 0; % all aggradation comes from feed

fIi(LogicArr_Aggr_Feed) = (1-PulseImpactRatio_Aggr(LogicArr_Aggr_Feed)) ...
    .* (alphaF_plusRandom.*NodeGSD.SurfPfi(LogicArr_Aggr_Feed) ...
    + (1 - alphaF_plusRandom)*NodeSedTrans.qbi(LogicArr_Aggr_Feed)) ...
    + PulseImpactRatio_Aggr(LogicArr_Aggr_Feed) .* NodeGSD.FeedPfi(LogicArr_Aggr_Feed) ;

% Transport Degrades and PulseFeed aggrades: NOT SURE  ABOUT THIS
% ONE RIGHT NOW>>>> Need to mix with subsurface?
% I think not, as feed goes all into transport
fIi(LogicArr_Degr_Feed) = pssi(LogicArr_Degr_Feed);

%         fIi(LogicArr_Degr_Feed) = (PulseImpactRatio_Aggr(LogicArr_Degr_Feed) + 1) .* pssi(LogicArr_Degr_Feed) ...
%             - (-PulseImpactRatio_Aggr(LogicArr_Degr_Feed) .* NodeGSD.FeedPfi(LogicArr_Degr_Feed));


%% Calculate transport rate per grain size class
for m = 1:gsd_MG
    sedTrans_qbxi(:,m) = NodeSedTrans.qb .* NodeSedTrans.qbi(:,m);
    feed_qbxi(:,m)     = NodeSedTrans.qb_feed .* NodeGSD.FeedPfi(:,m);
end

%% Calculating transport rate difference:

param_au_upwind   = RunParam.Upwind_au;
param_au_downwind = RunParam.Upwind_au - 1;
dqbxi = zeros(size(pssi));

switch RunParam.Upwind_order
    case 'first'
        % FIRST ORDER UPWIND SCHEME
        % https://en.wikipedia.org/wiki/Upwind_scheme\
        % first upstream node:
        %dqbxi(1,:) = (RunParam.au  .* sedTrans_qbxi(1,:)) - feed_qbxi(1,:)     ./ dx_array(1) ... % changed RunParam.au not to affect feed
        %    - (1 - RunParam.au) .* (sedTrans_qbxi(2,:) - sedTrans_qbxi(1,:))./ dx_array(1);    % changed to -(1-RunParam.au)
        dqbxi(1,:) = ((param_au_upwind   * sedTrans_qbxi(1,:)) - feed_qbxi(1,:))    ./ repmat(dx_array(1),1,gsd_MG)...
            + param_au_downwind * (sedTrans_qbxi(2,:) - sedTrans_qbxi(1,:))./ repmat(dx_array(1),1,gsd_MG);
        
        % all nodes in the middle:
        %dqbxi(2:end-1,:) = param_au_upwind   * (sedTrans_qbxi(2:end-1,:) - sedTrans_qbxi(1:end-2,:)) ./ repmat(dx_array(2:end-1),1,gsd_MG) ...
        %                 + param_au_downwind * (sedTrans_qbxi(  3:end,:) - sedTrans_qbxi(2:end-1,:)) ./ repmat(dx_array(2:end-1),1,gsd_MG); % changed to -(1-RunParam.au)
        % upstream and downstream:
        dqbxi(2:end-1,:) = param_au_upwind   * (sedTrans_qbxi(2:end-1,:)   - sedTrans_qbxi(1:end-2,:)) ./ repmat(dx_array(2:end-1),1,gsd_MG) ...
            + param_au_downwind * (sedTrans_qbxi(3:end,:) - sedTrans_qbxi(2:end-1,:)) ./ repmat(dx_array(2:end-1),1,gsd_MG);
        
        %dqbxi(end-1,:) = param_au_upwind  (sedTrans_qbxi(end-1,:) - sedTrans_qbxi(end-2,:)) ./ repmat(dx_array(end),1,gsd_MG);
        dqbxi(end,:) = param_au_upwind * (sedTrans_qbxi(end,:) - sedTrans_qbxi(end-1,:)) ./ repmat(dx_array(end),1,gsd_MG);
    case 'second'
        % SECOND ORDER UPWIND SCHEME https://en.wikipedia.org/wiki/Upwind_scheme
        % first node (use first order, rest second exept last 2 nodes):
        dqbxi(1,:) = ((param_au_upwind   * sedTrans_qbxi(1,:)) - feed_qbxi(1,:)    )./ repmat(dx_array(1),1,gsd_MG)...
            + param_au_downwind * (sedTrans_qbxi(2,:) - sedTrans_qbxi(1,:))./ repmat(dx_array(1),1,gsd_MG);
        
        dqbxi(2,:) = (param_au_upwind   * (  3 * sedTrans_qbxi(2,:) - 4 * sedTrans_qbxi(1,:)) + 1 * feed_qbxi(1,:)    )./(2 * repmat(dx_array(2),1,gsd_MG))...
            + param_au_downwind * (- 3 * sedTrans_qbxi(2,:) + 4 * sedTrans_qbxi(3,:) - 1 * sedTrans_qbxi(4,:))./(2 * repmat(dx_array(2),1,gsd_MG));
        
        dqbxi(3:end-2,:) = param_au_upwind   * (  3 * sedTrans_qbxi(3:end-2,:) - 4 * sedTrans_qbxi(2:end-3,:) + 1 * sedTrans_qbxi(1:end-4,:))./(2 * repmat(dx_array(3:end-2),1,gsd_MG))...
            + param_au_downwind * (- 3 * sedTrans_qbxi(3:end-2,:) + 4 * sedTrans_qbxi(4:end-1,:) - 1 * sedTrans_qbxi(  5:end,:))./(2 * repmat(dx_array(3:end-2),1,gsd_MG));
        
        % For this node first order:
        dqbxi(end-1,:) = param_au_upwind   * (sedTrans_qbxi(end-1,:) - sedTrans_qbxi(end-2,:)) ./ repmat(dx_array(end-1),1,gsd_MG)...
            + param_au_downwind * (sedTrans_qbxi(  end,:) - sedTrans_qbxi(end-1,:)) ./ repmat(dx_array(end-1),1,gsd_MG);
        
        dqbxi(end,:) = param_au_upwind * (sedTrans_qbxi(end,:) - sedTrans_qbxi(end-1,:)) ./ repmat(dx_array(end),1,gsd_MG);
    otherwise
        error('unknown upwinding scheme specified!')
end

% FEED IMPACT
dqbxi_feed = feed_qbxi ./ repmat(dx_array,1,gsd_MG);
dqbxi_feed(1,:) = zeros(1,gsd_MG); % node 1 feed is already in the upwinding scheme
dqbxi = dqbxi - dqbxi_feed;
sedTrans_array_dqbx_feed = NodeSedTrans.dqb-NodeSedTrans.dqb_feed;
sedTrans_array_dqbx_feed(1) = NodeSedTrans.dqb(1); % prevent feed impact on first node. Already in transport through upwind scheme

%keyboard

%% Calculating new ACTIVE LAYER GS frequencies:
for m = 1:gsd_MG
    NodeGSD.SurfPfi_new(:,m) = NodeGSD.SurfPfi(:,m) ...                               % BEFORE
        + 1./NodeSubs.ActiveLayerZ./(1 - RunParam.Lambda) ...
        .* (-dqbxi(:,m) + fIi(:,m) .* sedTrans_array_dqbx_feed(:)) ...
        .* RunParam.dt_actual.*RunParam.freqQ;    % + EXCHANGE IN?
    % ######## properly add pulse feed
    % OLD CARLES: psinew(i,k) = psi(i,k) + 1/La(i,1)/(1 - lambda)*(-dqbxi +  fIi(i,k)*dqbx(i,1))*dt*freqQ;
    if t_step > 1
        NodeGSD.SurfPfi_new(:,m) = NodeGSD.SurfPfi_new(:,m) ...
            - (1./[NodeSubs.ActiveLayerZ]) ...
            .* (NodeGSD.SurfPfi(:,m) - fIi(:,m)) ...
            .* ([NodeSubs.ActiveLayerZ] - [NodeSubs.ActiveLayerZ_old]);
        % OLD CARLES:psinew(i,k) = psinew(i,k) + (-psi(i,k) + fIi(i,k))*(1/La(i,1))*(La(i,1) - Laold(i,1));
    end
    assert(all(~isnan(NodeGSD.SurfPfi_new(:,m))),'non real Pfi!');
end
NodeGSD.SurfPfi_new(NodeGSD.SurfPfi_new<0) = deal(0);

%% Computing new Surface GSD parameters
SurfFsi = NaN(nodes_N,gsd_MG+1);
for n = 1:nodes_N
    % Recalc Fsi from Pfi
    SurfFsi(n,:) = [0, cumsum(NodeGSD.SurfPfi_new(n,:))];
    SurfFsi(n,:) = SurfFsi(n,:)./SurfFsi(n,end);
    
    assert(all(~isnan(SurfFsi(n,:))),'non real Fsi!');
end

% copy results to arrays that we return:
NodeSedTrans.fIi = fIi;
NodeSedTrans.dqbxi = dqbxi;

NodeGSD.SurfFsi = SurfFsi;


end

