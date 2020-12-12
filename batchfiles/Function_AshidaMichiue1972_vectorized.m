function [qbt,pbi] = Function_AshidaMichiue1972_vectorized(Dimean,pi,Dg,R,g,ustar)
%FUNCTION_ASHIDAMICHIUE1972 Bedload transport function as given by Ashida
%and Michiue 1972. Coded after the formulas provided by Garry Parker's ebook.
format long e
%keyboard

%% Allocating variables
gsd_MG = length(Dimean);

% per node and grain size class:
Wi_star     = zeros(1,gsd_MG);
tau_star_ratio = zeros(1,gsd_MG);
qbi_star    = zeros(1,gsd_MG);
qbi         = zeros(1,gsd_MG);
pbi         = zeros(1,gsd_MG);
Phi         = zeros(1,gsd_MG);
tau_i_star  = zeros(1,gsd_MG);
tau_ci_star = zeros(1,gsd_MG);
diff        = zeros(1,gsd_MG);

% dimensionless critical shear stress
tau_c_star = 0.05;

%% Shear stress:
%tau_sg_star = ustar(:).^2./(R*g.*Dg(:)./1000);
tau_i_star(:) = ustar(:).^2./(R*g.*Dimean(:)./1000);

%% hiding function:
% Modified version of Egiazaroff (1965) hiding function after Parker ebook:
% Di over Dsg:
diff(:) = Dimean(:)./Dg;

Phi(diff<= 0.4) = 0.843 .* diff(diff<= 0.4) .^ (-1); 
Phi(diff>0.4)   = (log(19) ./ (log(19 .* diff(diff>0.4)))) .^ 2;

% shear stress per grain size using hiding function:
tau_ci_star = Phi .* tau_c_star;

%% dimensionless bedload transport
tau_star_ratio(tau_i_star > tau_ci_star) = tau_ci_star(tau_i_star > tau_ci_star) ./ tau_i_star(tau_i_star > tau_ci_star);

Wi_star = 17 .* (1 - tau_star_ratio(tau_i_star > tau_ci_star)) .* (1 - sqrt(tau_star_ratio(tau_i_star > tau_ci_star)));
%qbi_star = 17 * ( tau_i_star - tau_ci_star) * (sqrt(tau_i_star)-sqrt(tau_ci_star));
% Sediment trasport rate of the ith fraction and total load (dimensional)
qbi(:) = (Wi_star(:).*pi(:).*ustar^3)./(R*g);
%qbi(:) = qbi_star(:) * pi(:) .* sqrt(R .* g .* Dimean(:)/1000) .* Dimean(:);
qbt = sum(qbi);

% Sediment transport rate frequencies
if qbt > 0
    pbi(:) = qbi(:)./qbt;
else
    pbi(:) = 0;
    qbt = 0;
end
