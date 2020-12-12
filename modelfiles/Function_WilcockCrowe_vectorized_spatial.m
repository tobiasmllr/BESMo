function [qbtN,pbiN,tausg_star,taussrgstar] = Function_WilcockCrowe_vectorized_spatial(Di,DimeanN,piN,DgN,...
    R,g,ustar,gsd_MG,nodes_N,t,tau_crit_factor)
% Compute total transport and bedload frequencies, using Wilcock and Crowe,
% 2003.
% The function needs six inputs to run:
% 1. Di: particle sizes. (in m??)
% 2. pi: bed surface frequencies.
% 3. Dg: geometrical mean size (m)
% 4. R: specifig gravity of sediment (Specific density???)
% 5. g: acceleration of gravity
% 6. ustar: shear velocity
format long e

%% SUBROUTINE: GSD Parameters: Dg and sg
% Grain classes
%MG = length(piN);

%% Allocating variables
%FacN         = zeros(nodes_N,MG+1);     
Wi_star     = zeros(nodes_N,gsd_MG);
pbiN        = zeros(nodes_N,gsd_MG);

psiN = log(transpose(repmat(Di,[1 nodes_N])))./log(2);

% Reconstructing cumulative function
FacN = [zeros(nodes_N,1), cumsum(piN,2)];

% sand content (%)
logicArray(:,:) = (psiN(:,1:end-1) <= 1) .* (psiN(:,2:end) > 1);
Fs = sum(...
    logicArray .* ((FacN(:,2:end) - FacN(:,1:end-1)) ...
    ./(psiN(:,2:end) - psiN(:,1:end-1)) ...
    .*(psiN(:,1:end-1) - 1) + FacN(:,1:end-1)...
    ),2);

zeroArray = zeros(1,nodes_N);
Fs(isnan(Fs)) = zeroArray(isnan(Fs));

%% dimensionless critical shear stress
taussrgstar = tau_crit_factor * (0.021 + 0.015*exp(-20*Fs));

%% Shear stresses
tausg_star = ustar(:).^2./(R*g.*DgN(:)./1000);

%% Hiding/Exposure coefficients
phisg0 = tausg_star./taussrgstar;

diff = DimeanN./repmat(DgN,[1 gsd_MG]);

b = 0.67./(1 + exp(1.5 - diff));

phi_i = repmat(phisg0,[1 gsd_MG]).*diff.^(-b);

Wi_star(phi_i < 1.35)   = 0.002*phi_i(phi_i < 1.35).^7.5;
Wi_star(phi_i >= 1.35)  = 14*(1 - 0.894./phi_i(phi_i >= 1.35).^0.5).^4.5;

% excess shearstress logical
excessTauLogic = tausg_star > taussrgstar;
% Sediment trasport rate of the ith fraction and total load (dimensional)
qbiN =  repmat(excessTauLogic,[1 gsd_MG]) .* (Wi_star .* piN .* repmat(ustar,[1 gsd_MG]).^3) ./ (R*g);

qbtN = sum(qbiN,2);

%qbtN MINIMUM 0
qbtN(qbtN<0)=0;

% Sediment transport rate frequencies
for n=1:nodes_N
    if qbtN(n) > 0
        pbiN(n,:) = qbiN(n,:)./qbtN(n);
    else
        pbiN(n,:) = 0;
        qbtN(n) = 0;
    end
end

%BREAK if problem
if isnan(sum(pbiN(:)))
    disp('isnan(pbiN)')
    disp(t)
    keyboard
end
%BREAK if problem
if isnan(sum(qbtN(:)))
    disp('isnan(qbtN)')
    disp(t)
    keyboard
end

%BREAK if problem
if isnan(sum(qbiN(:)))
    disp('isnan(qbi)')
    disp(t)
    keyboard
end

