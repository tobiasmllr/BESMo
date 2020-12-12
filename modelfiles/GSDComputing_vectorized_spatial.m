function [DgN,sgN,D90N,D50N,piN,Dimean,increaseSubtimestep] = GSDComputing_vectorized_spatial(FxiN,Dxi,gsd_MG,nodes_N,increaseSubtimestep)                                                       
format long e

%% Initial values of grain classes and frequencies.
%MG = size(FxiN,2)-1;

%% Allocating variables
psi_mi = zeros(nodes_N,gsd_MG);
piN = zeros(nodes_N,gsd_MG+1);

%for n=1:nodes_N
%    psi_i(n,:) = log(Dxi)/log(2);
%end
psi_i = log(transpose(repmat(Dxi(:),1,nodes_N)))/log(2);

%% Substracting frequencies smaller than 0
piN=FxiN(:,2:end)-FxiN(:,1:end-1);
piN(piN<0)=0;    % replace all values smaller than 0 with 0

%% Normalization of the GSD (making frequencies to sum 1)
% renormalized density function
sumpi = sum(piN,2);
for m=1:gsd_MG
    piN(:,m) = piN(:,m)./sumpi(:);
end

% Reconstructing cumulative function
FxiN = [zeros(nodes_N,1), cumsum(piN,2)];

%% Geometric mean diameter
% Mean size of each sediment class Dxmi, in mm
psi_mi(:,:) = (psi_i(:,1:gsd_MG) + psi_i(:,2:gsd_MG+1))./2;

Dimean = (2.^psi_mi);

psi_m = sum(piN.*psi_mi,2);

% Geometric mean diameter, in mm
DgN = 2.^psi_m;

if any(isnan(DgN))
    warning('non real Dg!');
    increaseSubtimestep = true;
end

for i = 1:nodes_N
    psi_mDiff(i,:) = psi_m(i) - psi_mi(i,:);
end

%% Geometric standard deviation, (-)
sigma2 = sum(piN.*psi_mDiff.^2,2);

% Geometric standard deviation, (-)
sgN = 2.^sqrt(sigma2);

%% Calculation of D90, in mm. Linear interpolation in psi variables
logicArray          = (FxiN(:,1:end-1) <= 0.9) .* (FxiN(:,2:end) > 0.9);
logicGreaterZero    = FxiN(:,:) <= 0;
psi90 = nansum(logicArray .* ...
    (...
    psi_i(:,1:end-1) + ...
    (psi_i(:,2:end)   - psi_i(:,1:end-1)) ./ ...
    (FxiN(:,2:end)    - FxiN(:,1:end-1)) .* ...
    (0.9              - FxiN(:,1:end-1))...
    ),2);
%posi = find (logicArray == 2);
%psi90 = psi_i(posi)+ (psi_i(posi+1) - psi_i(posi))/(Fxi2(posi+1) - Fxi2(posi))*(0.9 - Fxi2(posi));
D90N = 2.^psi90;

%% Calculation of D50, in mm. Linear interpolation in psi variables
logicArray          = (FxiN(:,1:end-1) <= 0.5) .* (FxiN(:,2:end) > 0.5);
logicGreaterZero    = FxiN(:,:) <= 0;
psi50 = nansum(logicArray .* ...
    (...
    psi_i(:,1:end-1) + ...
    (psi_i(:,2:end)   - psi_i(:,1:end-1)) ./ ...
    (FxiN(:,2:end)    - FxiN(:,1:end-1)) .* ...
    (0.5              - FxiN(:,1:end-1))...
    ),2);
D50N = 2.^psi50;
