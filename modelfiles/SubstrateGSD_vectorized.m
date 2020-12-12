function [pssi_new,NodeSubs] = SubstrateGSD_vectorized(...
    RunParam,NodeGSD,NodeGeom,NodeSubs,...
    fIi,LstrMat,gsd_MG)

% [pssi_new,Msj_new,eta_st_new,pssi_st_new] = SubstrateGSD_vectorized(...
%     MG,nLa,D90s_new,D90s,eta_new,eta,Msj,deta_s,LstrMat,eta_st,zbase(n),...
%     pssi_st,fIi,psinew,gsd_feedPfi,PulseImpactRatio)
%Debugging Tobias
% MG = gsd_MG
%NodeGSD.SurfPfi,NodeGSD.SurfDg,NodeGSD.SurfSDg,NodeGSD.SurfD90,NodeGSD.SurfDimean
%NodeGeom.Elev_new

format long g

%%
psinew  = NodeGSD.SurfPfi_new;
pssi_st = NodeSubs.Pfi;
eta_st  = NodeSubs.Elev;
eta     = NodeGeom.Elev;
eta_new = NodeGeom.Elev_new;
deta_s = RunParam.deta_s;
nodesN = size(eta,1);

% Msj is the first surface layer position (below active layer).
% In new geometry, the active layer is always at index 1 (Msj-1), so the first subsurface
% layer is always at index 2 (Msj)
Msj = ones(nodesN,1) * 2;
Msj_new = ones(nodesN,1) * 2;

D90s    = NodeGSD.SurfD90;
D90s_new= NodeGSD.SurfD90_new;
nLa     = RunParam.nLa;

%zbase is not useful in this new geometry. replaced with checks of the
%depth of first subsurface layer below interface
%zbase = NodeGeom.Elev-NodeSubs.ActiveLayerZ-RunParam.SubsurfaceDepth;

% How the older version stores subsurface data per node:
% 19 ...
% 18  ~~~~~~ Water / Air (ignored)
% 17  ~~~~~~ Water / Air (ignored)
% 16  ~~~~~~ Msj+2: Water / Air (ignored)           | ERROR in original logic? SHOULD CONSIDER zbase(n)?
% 15  ###### Msj+1: Active Layer (transport here)   | Elevation above base:  s
% 14  oooooo Msj+0: Subsurface Layer 1 (surface)    | top: eta_inter_new :
% 13  xxxxxx Msj-1: Subsurface Layer 2              | Elevation above base:
% 12  xxxxxx Msj-2: Subsurface Layer 3 etc.         | Elevation above base:
% 11 ...

% New version:
% -1 Temporarily possible to move active layer here, but then we circshift back to index 0
% 0 Same circshift back to index 0
% 1  ###### Msj-1: Active Layer (transport here)    | Elevation above base:
% 2  oooooo Msj+0: Subsurface Layer 1               | eta_inter_new       :
% 3  xxxxxx Msj+1: Subsurface Layer 2               | Elevation above base:
% 4  xxxxxx Msj+2: Subsurface Layer 3 etc.          | Elevation above base:
% 5 ...

%% Updating the stratigraphy for the new time-step.

La_str_new = nLa*D90s_new./1000;
eta_inter_new = eta_new - La_str_new;

La_str = nLa*D90s./1000;
eta_inter = eta-La_str;
%diff_eta = eta_inter_new - eta_inter;

% unrealistic check now, replaced by looking for NaN later!:
%assert(all(eta_inter_new > eta_st(:,end)),'Subsurface elevation not deep enough! Maybe increase param_alfa_LZd')

%% if eta_inter_new < eta_inter
% FIRST we check if the new interface between active layer and first
% subsurface layer is moving
%listDegradingNodes = find(eta_inter_new < eta_inter);
listDegradingNodes = find(eta_inter_new < eta_st(:,2));
if any(listDegradingNodes);
    % Found some nodes, where the new interface is within the subsurface!
    for nIter=1:length(listDegradingNodes)
        % loop only through nodes that are changing, so that we can
        % break out of subsurface loop and save some computations,
        % as most times the active layer will still be at the surface.
        % No need to iterate through everything!
        n = listDegradingNodes(nIter);
        
        lyr = 2;
        while lyr < LstrMat
            % Check if we ran out of subsurface:
            assert(~isnan(eta_st(n,lyr+1)),'Error: found NaN value in subsurface stratigraphy that affects active layer');
            
            if (eta_inter_new(n) <= eta_st(n,lyr)) && (eta_inter_new(n) > eta_st(n,lyr+1))
                % set new surface layer index to layer above to where this is true
                Msj_new(n) = lyr;
                % lyr   is new first subsurface
                % lyr-1 is new surface layer
                % lyr-2 is new active layer
                
                break % stop iteration
            else
                lyr = lyr + 1;
            end
        end
    end
end


%% elseif eta_inter_new > eta_inter
listAggradingNodes = find(eta_inter_new > eta_st(:,2));
if any(listAggradingNodes);
    % at all nodes where the first subsurface layer was moved due to aggradation
    % we have to find the new activelayer position:
    
    for nIter=1:length(listAggradingNodes)
        n = listAggradingNodes(nIter);
        %beta = eta_inter_new(n) - zbase(n) - RunParam.SubsurfaceDepth;
        
        % trying to get rid of zbase:
        beta = eta_inter_new(n) - eta_st(n,3);
        %keyboard
        
        %L = eta_inter_new - zbase(n) - ((Msj - 2)*deta_s);
        Msj_new(n) = Msj(n) - floor(beta/deta_s);
        %Msj_new = Msj + floor(L/deta_s);
    end
end
% else
%     Msj_new = Msj;
% end


%% Setting the dimensions of the stratigraphy
% take over old values and overwrite on demand:
eta_st_new  = eta_st;
pssi_st_new = pssi_st;
pssi_new(:,:) = pssi_st(:,2,:);

%pssi_temp = zeros(1,gsd_MG);
%logicCounter = [1:1:LstrMat];

%% Array geometry
% circshift array, so that surface layer is index 1, second subsurface 2
% etc. The new position of the active layer Msj_new should be 0 everywhere!
shift_layerDist = Msj-Msj_new;
listShiftNodes = find(shift_layerDist);
if any(listShiftNodes)
    for nIter=1:length(listShiftNodes)
        n = listShiftNodes(nIter);
        % move the elements to new positions in array
        eta_st_new(n,:) = circshift(eta_st_new(n,:),shift_layerDist(n),2);
        pssi_st_new(n,:,:)= circshift(pssi_st_new(n,:,:),shift_layerDist(n),2);
        
        if shift_layerDist(n) < 0
            % if the surface is degragading:
            bottomLyr = LstrMat + shift_layerDist(n)-1;
            
            % overwrite all elements that are now at the bottom with NaN
            eta_st_new(n,bottomLyr:end) = deal(NaN());
            pssi_st_new(n,bottomLyr:end,:) = deal(NaN());
        else
            % if the surface has accumulated:
            topLyr = shift_layerDist(n) + 2;
            
            % overwrite all elements that are now at the top with NaN, we
            % will find the proper values later
            eta_st_new(n,1:topLyr) = deal(NaN());
            pssi_st_new(n,1:topLyr,:) = deal(NaN());
        end
    end
end

%%
%keyboard
%% 1. Aggradation Case
% trying to get rid of zbase:
%beta_old    = eta      - eta_st(:,2);
%beta_new    = eta_new  - eta_st(:,2);

%if eta_inter_new > eta_st(:,1)
if any(listAggradingNodes);
    for nIter=1:length(listAggradingNodes)
        n = listAggradingNodes(nIter);
        
        beta_old = eta_inter(n) - eta_st(n,3);
        beta_new = eta_inter_new(n) - eta_st(n,3);
        
        % In all aggrading cases:
        % new Active Layer
        eta_st_new(n,1) = eta_new(n);       % Elevation
        
        % get frequencies from transport
        pssi_st_new(n,1,:) = psinew(n,:);   % Frequencies (pulse added earlier):
        
        if Msj_new(n) == Msj(n)
            %% 1.1 Aggradation is not high enough to add another computational point in the stratigraphy
            
            % Uppermost substrate Layer
            % Elevation
            eta_st_new(n,2)     = eta_new(n) - La_str_new(n);
            
            % Elevation difference:
            % beta is the difference of the surface to the first subsurface
            % layer (so substracting active layer thickness La_str, and
            % Interface elevation (old: (Msj -2)*detas
            %             beta_old(n)    = eta(n)      -La_str(n)     - (Msj(n)     - 2)*deta_s - zbase(n);
            %             beta_new(n)    = eta_new(n) - La_str_new(n) - (Msj_new(n) - 2)*deta_s - zbase(n);
            %beta_old(n)    = eta(n)      -La_str(n)     - zbase(n) - RunParam.SubsurfaceDepth;
            %beta_new(n)    = eta_new(n) - La_str_new(n) - zbase(n) - RunParam.SubsurfaceDepth;
            
            
            %keyboard
            
            % Similar stuff from earlier, I'm still not sure what beta
            % represents. Distance of interface from a layer boundary?
            %beta = eta_inter_new(n) - zbase(n) - RunParam.SubsurfaceDepth;
            %L = eta_inter_new - zbase(n) - ((Msj - 2)*deta_s);
            %Msj_new(n) = Msj(n) - floor(beta/deta_s);
            
            % Frequencies of the layer underneath the Active Layer
            pssi_new(n,:) = 1/beta_new.*(fIi(n,:).*(beta_new - beta_old) + transpose(squeeze(pssi_st(n,3,:) .* beta_old)));
            
            % Store new frequencies in first substrate layer
            pssi_st_new(n, 2,:)  = pssi_new(n,:);
            
        else
            %% 1.2 Aggradation is high enough to add another computational point in the stratigraphy
            
            % lowest new layer:
            lyrBot = Msj(n) + shift_layerDist(n);
            
            % fill up subsurface that is not changing
            eta_st_new(n,lyrBot+1:end) = eta_st(n,3:end-shift_layerDist(n));
            
            % get list of layers that have to change
            % except active layer (index 1), except bottom aggrading layer (index lyrBot)
            % and except first surface layer (index 2):
            lyrArray = [lyrBot-1:-1:Msj(n)+1];
            
            if ~isempty(lyrArray)
                % if there is more than one layer added (more than the lyrBot and first substrate lyr)
                
                % 1.2.2. More than one point is added in the stratigraphy.
                % These are the central layers between interface old
                % subsurf/newsubsurf (handled in 1.2.1)
                % and active layer/topsubsurf (handled in 1.2.1)
                
                for lyr=1:length(lyrArray)
                    % lyr is an index that signifies how many layer thicknesses
                    % higher than the old subsurface this layer is. The lyr=1
                    % is the case wiyh lyrBot below...
                    
                    % get new elevation from elevation of old first
                    % subsurface layer (always lyr index 3) with distance
                    % in lyr*deta_s steps
                    eta_st_new(n,lyrArray(lyr))     = eta_st(n,Msj(n)+1) + deta_s * (lyr + 1);
                    pssi_st_new(n,lyrArray(lyr),:)  = fIi(n,:);
                end
            end
            
            %% 1.2.1. Top and bottom addition:
            %% Lowest point to subsurface!
            %                 % Substrate layer -------------------------------------
            %                 eta_st_new(logicCounter < Msj_new-1)      = eta_st(logicCounter < Msj_new-1);
            %                 pssi_st_new(logicCounter < Msj_new-1,:)   = pssi_st(logicCounter < Msj_new-1,:);
            %                 beta_old(n)                                    = eta - La_str - (Msj - 2)*deta_s - zbase(n);
            
            eta_st_new(n,lyrBot)     = eta_st(n,Msj(n)+1)+ deta_s;
            
            % get frequencies from layer below this!
            pssi_new(n,:)               = 1/deta_s*(fIi(n,:).*(deta_s - beta_old) + transpose(squeeze(pssi_st(n,Msj(n)+1,:)*beta_old)));
            pssi_st_new(n,lyrBot,:)  = pssi_new(n,:);
            
            %% New Uppermost substrate layer (always index 2 in new geom)
            eta_st_new(n,2)       = eta_new(n)-La_str_new(n);
            
            % get frequencies from bed surface exchange
            pssi_st_new(n,2,:)    = fIi(n,:);
            pssi_new(n,:) = pssi_st_new(n,2,:);
            
            % test if there is a layer that's too high
            assert(all(min(diff(eta_st_new(n,2:lyrBot+1)))+deta_s+10E-6>0),'some subsurface layers are higher than allowed!')
        end
    end
end
%else
%% 2. Degradation Case
if any(listDegradingNodes);
    for nIter=1:length(listDegradingNodes)
        n = listDegradingNodes(nIter);
        %         % Uppermost substrate Layer
        %         eta_st_new(logicCounter == Msj_new)         = eta_new - La_str_new;
        %         pssi_st_new(logicCounter == Msj_new,:)      = pssi_st(logicCounter == Msj_new,:);
        %         pssi_new(1,:)                               = pssi_st(logicCounter == Msj_new,:);
        
        % Uppermost substrate Layer
        eta_st_new(n,2)         = eta_new(n) - La_str_new(n);
        pssi_st_new(n,2,:)      = pssi_st(n,2,:);
        pssi_new(n,:)           = pssi_st(n,2,:);
        
        % Active Layer
        eta_st_new(n,1)     = eta_new(n);
        pssi_st_new(n,1,:)  = psinew(n,:);
    end
end


% I hope this is ok to do for the whole array / Tobias:
%pssi_st_new(pssi_st_new < 0) = 0;

%% overwrite old values with new results
NodeSubs.Pfi  = pssi_st_new;
NodeSubs.Elev = eta_st_new;


% for n=1:nodesN
%     differences_new=diff(eta_st_new(n,1:10),1,2);
%     assert(all(differences_new(:)+deta_s+10E-6>0),'some subsurface layers are higher than allowed!')
%     keyboard
% end

%% Debugging
% if(isnan(sum(eta_st_new(1:Msj_new + 1))))
%     error('Error: isnan(sum(eta_st_new(1:Msj_new + 1)))')
% end
