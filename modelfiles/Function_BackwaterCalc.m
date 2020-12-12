function [ustar_array,waterDepth_final,froude_array_final] = Function_BackwaterCalc(qw,...
    slope_array, BedElev, gsd_arraySurfD90_new,...
    dx_array,nodes_loc,param_nk,param_alfa_r,...
    const_rhoW,const_g,backwaterparam)
format long e

%% fake input to understand the script without calling it from the model:
debug = false;
if debug
    qw                      = 0.0650;
    % flat downslope:
    slope_array             = ones(41,1) .* 0.05;
    slope_array(10:18)      = -0.01;
    
    % bumpy:
%     slope_array             = transpose([0.01 0.01 0.01 0.01 0.01 0.05 0.05 0.05 ...
%         0.05 0.05 0.07 0.07 0.07 0.07 0.05 0.05 ...
%         0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 ...
%         0.03 0.03 0.03 0.03 0.03 0.03 0.07 0.07 ...
%         0.01 0.01 0.01 0.01 0.001 0.001 0.001 0.001 ...
%         0.001]);
    
    gsd_arraySurfD90_new    = ones(41,1) .* 23.425;
    dx_array                = ones(41,1) .* 0.5;
    nodes_loc               = zeros(length(dx_array),1);
    nodes_loc(2:end)        = cumsum(dx_array(1:end-1));
    param_nk                = 2;
    param_alfa_r            = 8.1;
    const_rhoW              = 1000;
    const_g                 = 9.81;
    
    BedElev                 = NaN(length(dx_array),1);
    BedElev(end) = 0; % outflow
    for n = length(BedElev)-1:-1:1
        BedElev(n) = BedElev(n+1) + tan(slope_array(n))*dx_array(n);
    end
    
    backwaterparam.iter_max         = 1000;     % maximum number of iterations to reach iter_diffgood:
    backwaterparam.iter_diffgood    = 0.00001;   % stop iterating waterdepth approximation if difference between calculations is smaller than this value
    backwaterparam.waterdepth_min   = 0.00001;   % calc finer spatial grid if waterdepth is calculated to be smaller than value. maybe calculate this threshold realistically?
    backwaterparam.froudefact_min   = 0.8;      % calc finer spatial grid if upstream decrease of froudenumber in subcritical condition is calculated to be smaller than this value multiplied by downstream Froude number.
    backwaterparam.subGridScalemax  = 65536;    % maximum interpolation factor for calculating subgrids. Error if tripped
    backwaterparam.Fr_threshold     = 0.9;      % Froude Number Threshold for using supercritical:
end

%% Backwater curve calcutation
% Nikuradse coefficient, in m
ks = param_nk.*gsd_arraySurfD90_new(:)./1000;

%% Set up arrays
nodes_N = length(slope_array);
waterDepth_normal   = zeros(nodes_N,1);
taub                = NaN(nodes_N,1);
waterDepth_final    = NaN(nodes_N,1);
froude_array_final  = NaN(nodes_N,1);

% Critical water depth, in m
waterDepth_critical = (qw^2/const_g)^(1/3);

% Normal Water depth Calculation (Manning-strickler), in m
waterDepth_normal=((ks(end)^(1/3)*qw^2)/(param_alfa_r^2*const_g*slope_array(end)))^(3/10);

% We start calculation at the outlet (downstream boundary condition)
waterDepth_downstreamBoundary = waterDepth_normal(end);

% Froude Number Threshold for using supercritical:
%Fr_threshold = 0.9;

% Prepare parameters for calculation:
doLoop = true;
subGridScale=1;
subGridFactor=4;

% for first try use original resolution:
slope_array_sample = slope_array;
ks_sample = ks;
dx_array_sample = dx_array;
nodes_loc_original = nodes_loc;

while doLoop
    [froude_array_sample,frictionSlope_array_sample, waterDepth_final_sample, redoWithSubStep] = BackwaterArray(qw,slope_array_sample,ks_sample,...
        dx_array_sample,param_alfa_r,const_g,...
        waterDepth_downstreamBoundary,waterDepth_critical,backwaterparam);
    if redoWithSubStep
        % There was a problem with resolution and we have to increase the
        % resolution for the backwater calculation
        subGridScale       = subGridScale * subGridFactor;
        if subGridScale > backwaterparam.subGridScalemax
            warning('maximum subGridScale tripped. Something is going wrong, is backwaterparam.froudenumb_min or backwaterparam.waterdepth_min too small?');
            %keyboard
            %figure;plot(waterDepth_final_sample);hold on; plot(froude_array_sample)
        end
        grid_sample = interp1(1:nodes_N,nodes_loc_original,1:(1/subGridScale):length(nodes_loc_original));
        
        %grid_sample         = 1:(dx_array/subGridScale):nodes_N;
        slope_array_sample  = interp1(nodes_loc_original, slope_array,grid_sample);
        ks_sample           = interp1(nodes_loc_original, ks,grid_sample);
        dx_array_sample               = NaN(length(grid_sample),1);
        dx_array_sample(1:end-1)      = grid_sample(2:end) - grid_sample(1:end-1);
        dx_array_sample(end)          = grid_sample(end -1);
    else
        %figure;plot(waterDepth_final_sample);
        % We don't have to recalculate, but might have to downscale the results:
        if subGridScale > 1
            [~,indices] = ismember(nodes_loc_original,grid_sample);
            %waterDepth_final     = waterDepth_final_sample(indices);
            taub_sample         = const_rhoW .* const_g .* waterDepth_final_sample .* frictionSlope_array_sample; % Shear stress
            
            %% Get mean for supersampled values
            for i=1:length(indices)
                taub(i)                 = taub_sample(indices(i));
                waterDepth_final(i)     = waterDepth_final_sample(indices(i));
                froude_array_final(i)   = froude_array_sample(indices(i));
            end
        else
            waterDepth_final        = waterDepth_final_sample;
            froude_array_final      = froude_array_sample;
            taub = const_rhoW .* const_g .* waterDepth_final .* frictionSlope_array_sample; % Shear stress
            
        end
        break
    end
end
%%
if ~isreal(waterDepth_final)
    error('waterDepth_final gave non-real number!')
end

if debug
    %% original-node profile:
    figure;plot(BedElev,'-*k');hold on; plot(BedElev+waterDepth_final,'-*b'); hold on; plot(froude_array_final,'-*')
    
    %% supersample profile:
    BedElev_sample   = NaN(length(dx_array_sample),1);
    BedElev_sample(end) = 0; % outflow
    for n = length(BedElev_sample)-1:-1:1
        BedElev_sample(n) = BedElev_sample(n+1) + tan(slope_array_sample(n))*dx_array_sample(n);
    end
    figure;plot(BedElev_sample,'-k');hold on; plot(BedElev_sample+waterDepth_final_sample,'-b'); %hold on; plot(froude_array_sample,'-')
end

%% Calculate resulting shear stress
% set to zero in nodes with negative bedslope:
ustar_array = (slope_array > 0) .* sqrt(taub./const_rhoW); % Shear velocity

% debugging:
% if any(slope_array < 0)
%     keyboard
% end

% don't set to zero in nodes with negative bedslope:
%ustar_array = sqrt(taub./const_rhoW); % Shear velocity
end

function [froudeArray,frictionSlopeHere,waterDepth, redoWithSubStep] = BackwaterArray(qw,slope_array,ks,...
    dx_array,param_alfa_r,const_g,...
    waterDepth_downstreamBoundary,waterDepth_critical,backwaterparam)
% arrays
nodes_N             = length(slope_array);
frictionSlopeHere   = NaN(nodes_N,1);
waterDepth          = zeros(nodes_N,1);
froudeArray         = zeros(nodes_N,1);

stopLoop1 = false;
redoWithSubStep = false;
while ~stopLoop1
    % DOWNSTREAM END
    if waterDepth_downstreamBoundary >= waterDepth_critical
        waterDepth(end) = waterDepth_downstreamBoundary;
    else
        waterDepth(end) = waterDepth_critical + 0.001; % carles has it a bit higher +0.1;
    end
    frictionSlopeHere(end) = (ks(end)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth(end)^(10/3));
    froudeArray(end) = qw/const_g^(1/2)/waterDepth(end)^(3/2);
     
    for n = nodes_N-1:-1:1
        %% iterate through nodes backwards, from downstream to upstream:
        frictionSlope_DS = frictionSlopeHere(n+1);
        Slope_DS = slope_array(n+1);
        waterDepth_DS = waterDepth(n+1);
        
        % Calculate initial approximation of Froude number by using
        % downstream waterdepth / boundary waterdepth for downstream end:
        Fr1 = qw/const_g^(1/2)/waterDepth_DS^(3/2);
        
        if Fr1 > backwaterparam.Fr_threshold
            %% SUPERCRITICAL flow: normal flow approximation
            waterDepth(n) = ((ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*slope_array(n)))^(3/10);
            frictionSlopeHere(n) = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth(n)^(10/3));
            froudeArray(n) = qw/const_g^(1/2)./waterDepth(n).^(3/2);
        else
            %% SUBCRITICAL flow
            fSlope_i = frictionSlope_DS;
            Froude_i = Fr1;
            stopLoop2 = false;
            iter = 1;
            
            while ~stopLoop2
                dHdX_i = (Slope_DS - fSlope_i)/(1 - Froude_i^2);
                waterDepth_i = waterDepth_DS - dHdX_i*dx_array(n);
                fSlope_i = (ks(n)^(1/3)*qw^2)/(param_alfa_r^2*const_g*waterDepth_i^(10/3));
                Froude_i = qw/const_g^(1/2)/waterDepth_i^(3/2);
                
                % check how much the additional iteration changed the
                % value, if below threshold, then keep:
                if iter>1 && abs(dHdX_i - dHdX_old) < backwaterparam.iter_diffgood
                    % found value, go to next node
                    stopLoop2  = true;
                    break
                end
                
                % check if something went out of bound:
                if iter > backwaterparam.iter_max || ...
                        ~isreal(waterDepth_i * dHdX_i * fSlope_i * Froude_i)
                    % Problem A) too many iterations to find a value within iter_diffgood between two successive iterations
                    % Problem B) any of waterDepth_i, dHdX_i, fSlope_i, Froude_i is nonreal/imaginary number (happens if negative number^(10/3))
                    
                    % redo with higher resolution
                    stopLoop2 = true;
                    stopLoop1 = true;
                    redoWithSubStep = true;
                    break
                end
                
                % iterate
                iter = iter + 1;
                dHdX_old = dHdX_i;
            end
            
            % keep values:
            waterDepth(n) = waterDepth_i;
            frictionSlopeHere(n) = fSlope_i;
            froudeArray(n) = Froude_i;
        end
        
        if redoWithSubStep || ...
                waterDepth(n) < backwaterparam.waterdepth_min || ...         % OR Problem C: waterDepth_i is very small / negative
                n ~= nodes_N && (froudeArray(n) < (froudeArray(n+1) * backwaterparam.froudefact_min))            % OR Problem D: abrupt decrease of Froude number compared to downstream
            
            % Problem C) waterDepth_i is very small / negative
            % FROUDE CONDITION:
            % Allowed: smooth upstream decrease of Froude: supercritical --> critical --> subcritical
            % Allowed: abrupt upstream increase of Froude: subcritical --> hydr.jump --> supercritical
            % Problem D) abrupt upstream decrease of Froude: supercritical --> subcritical
            
            % redo with higher resolution
            stopLoop1 = true;
            redoWithSubStep = true;
            break % spatial loop
        end
    end
    break
end

end